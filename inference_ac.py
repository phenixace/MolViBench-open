#!/usr/bin/env python
"""
MolViBench - Agent Collaboration (AC)
=======================================

Two-agent framework (Coder + Tester):
    A1 (Coder): Generates code from the task description
    A2 (Tester): Designs test plan, executes code, judges chemical correctness,
                 and provides targeted feedback if it fails

Usage:
    python inference_ac.py --base-url https://api.openai.com/v1 --api-key sk-xxx --model gpt-4o
    python inference_ac.py ... --model-tester gpt-4o-mini --max-rounds 3 --lang cn --resume
"""

import argparse
import asyncio
import csv
import json
import os
import sys
import re
import time
import tempfile
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from evaluate.molecules import SINGLE_MOLECULES


# ============================================================
# Code Execution (reused from react agent)
# ============================================================

def execute_code(code: str, test_smiles: str, timeout: int = 30) -> str:
    """Execute code in subprocess, return output value or error traceback."""
    code_dir = tempfile.mkdtemp()
    code_path = os.path.join(code_dir, "user_code.py")
    runner_path = os.path.join(code_dir, "runner.py")

    with open(code_path, "w", encoding="utf-8") as f:
        f.write(code)

    runner = f'''import sys, json, traceback, importlib.util

spec = importlib.util.spec_from_file_location("user_code", {repr(code_path)})
mod = importlib.util.module_from_spec(spec)
try:
    spec.loader.exec_module(mod)
except Exception as e:
    tb = traceback.format_exc()
    print("__RESULT__" + json.dumps({{"ok": False, "error": f"Import error: {{type(e).__name__}}: {{e}}\\n{{tb}}"}}, default=str))
    sys.exit(0)

if not hasattr(mod, "level_function"):
    print("__RESULT__" + json.dumps({{"ok": False, "error": "Function level_function not defined."}}))
    sys.exit(0)

try:
    output = mod.level_function({repr(test_smiles)})
    print("__RESULT__" + json.dumps({{"ok": True, "output": repr(output), "type": type(output).__name__}}, default=str))
except Exception as e:
    tb = traceback.format_exc()
    print("__RESULT__" + json.dumps({{"ok": False, "error": f"{{type(e).__name__}}: {{e}}\\n{{tb}}"}}, default=str))
'''

    with open(runner_path, "w", encoding="utf-8") as f:
        f.write(runner)

    try:
        proc = subprocess.run(
            [sys.executable, runner_path],
            capture_output=True, text=True, timeout=timeout,
            env={**os.environ, "PYTHONIOENCODING": "utf-8"},
        )
        for line in proc.stdout.split("\n"):
            if line.startswith("__RESULT__"):
                data = json.loads(line[len("__RESULT__"):])
                if data["ok"]:
                    return f"SUCCESS | Output: {data['output']} | Type: {data['type']}"
                else:
                    return f"ERROR | {data['error']}"
        return f"NO_RESULT | stdout: {proc.stdout[:500]} | stderr: {proc.stderr[:500]}"
    except subprocess.TimeoutExpired:
        return f"TIMEOUT | Execution timed out ({timeout}s)."
    except Exception as e:
        return f"EXEC_ERROR | {e}"
    finally:
        import shutil
        shutil.rmtree(code_dir, ignore_errors=True)


# ============================================================
# Code Extraction
# ============================================================

def extract_code(text: str) -> str:
    """Extract Python code from LLM response."""
    if not text:
        return ""
    pattern = r"```(?:python)?\s*\n(.*?)```"
    matches = re.findall(pattern, text, re.DOTALL)
    if matches:
        for m in matches:
            if "level_function" in m:
                return m.strip()
        return matches[-1].strip()
    if "def level_function" in text:
        lines = text.split("\n")
        code_lines = []
        in_code = False
        for line in lines:
            if line.startswith(("import ", "from ", "def ", "#")) or in_code:
                in_code = True
                code_lines.append(line)
        if code_lines:
            return "\n".join(code_lines).strip()
    return text.strip()


# ============================================================
# API Call Helper
# ============================================================

async def call_llm(
    session, base_url: str, api_key: str, model: str,
    messages: list, temperature: float = 0.0,
    max_tokens: int = 4096, timeout: int = 120,
) -> dict:
    """Call OpenAI-compatible chat API. Returns {"content": ...} or {"error": ...}."""
    import aiohttp
    import ssl as _ssl

    url = f"{base_url.rstrip('/')}/chat/completions"
    headers = {"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"}
    payload = {
        "model": model,
        "messages": messages,
        "temperature": temperature,
        "max_tokens": max_tokens,
    }

    ssl_ctx = _ssl.create_default_context()
    ssl_ctx.check_hostname = False
    ssl_ctx.verify_mode = _ssl.CERT_NONE

    try:
        async with session.post(
            url, json=payload, headers=headers,
            timeout=aiohttp.ClientTimeout(total=timeout),
            ssl=ssl_ctx,
        ) as resp:
            if resp.status != 200:
                text = await resp.text()
                return {"error": f"HTTP {resp.status}: {text[:300]}"}
            data = await resp.json()
        return {"content": data["choices"][0]["message"].get("content", "")}
    except asyncio.TimeoutError:
        return {"error": f"Timeout ({timeout}s)"}
    except Exception as e:
        return {"error": str(e)}


# ============================================================
# Prompts
# ============================================================

# -- A1 (Coder) prompts --

A1_SYSTEM = """You are an expert cheminformatics programmer. You write Python functions using RDKit to solve molecular science tasks.

Rules:
- The function MUST be named `level_function`
- Include all necessary imports inside or above the function
- Handle invalid inputs gracefully (return None)
- Only use: rdkit, numpy, pandas, scikit-learn, matplotlib, selfies
- Return ONLY the Python code, no explanations"""

A1_GENERATE_USER = """Write a Python function `level_function` to solve this task:

Task: {question}

Requirements:
1. Function name must be `level_function`
2. Use the RDKit library
3. Include error handling (return None for invalid input)
4. Provide the complete Python code."""

A1_FIX_USER = """Your previous code for this task was tested and FAILED.

Task: {question}

Your previous code:
```python
{code}
```

Test results:
{exec_results}

Tester feedback:
{feedback}

Fix the code. IMPORTANT:
- Focus on fixing the SPECIFIC issue reported by the tester
- Do NOT change the overall chemical logic unless the tester explicitly says it is wrong
- Return the complete fixed Python code"""

# -- A2 (Tester) prompts --

A2_SYSTEM = """You are a chemistry QA expert who validates cheminformatics code. You evaluate whether code outputs are chemically correct and reasonable.

You must respond ONLY with valid JSON. No other text."""

A2_DESIGN_USER = """Given a molecular computing task, design a test plan BEFORE seeing any code.

Task: {question}

Available test molecules (SMILES):
{molecules}

Design your test plan as JSON with these fields:
{{
  "test_smiles": ["smiles1", "smiles2"],
  "expected_type": "int|float|str|list|dict|bool|tuple",
  "value_description": "what a correct output should look like",
  "reasonable_range": "e.g. 100-1000 for MW, true/false for boolean, etc.",
  "failure_modes": ["list of things that could go wrong"]
}}

Respond with ONLY the JSON object."""

A2_VERIFY_USER = """You designed this test plan for the task:

Task: {question}

Test plan:
{test_plan}

The code was executed with these results:
{exec_results}

Evaluate:
1. Does the code run without errors on all test molecules?
2. Is the output type correct (expected: {expected_type})?
3. Is the output chemically reasonable? {value_description}
4. Are there any obvious chemical logic errors?

Respond with ONLY a JSON object:
- If PASS: {{"verdict": "PASS"}}
- If FAIL: {{"verdict": "FAIL", "reason": "specific description of what is wrong and how to fix it"}}"""


# ============================================================
# Question Loading
# ============================================================

def load_questions(csv_path: str) -> list:
    questions = []
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            questions.append(row["question"])
    return questions


# ============================================================
# Helper: parse JSON from LLM response (tolerant)
# ============================================================

def parse_json_response(text: str) -> dict | None:
    """Try to extract a JSON object from LLM text."""
    if not text:
        return None
    # Try direct parse
    text = text.strip()
    # Strip markdown code fences if present
    if text.startswith("```"):
        text = re.sub(r"^```(?:json)?\s*\n?", "", text)
        text = re.sub(r"\n?```\s*$", "", text)
        text = text.strip()
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass
    # Try to find JSON object in text
    match = re.search(r"\{.*\}", text, re.DOTALL)
    if match:
        try:
            return json.loads(match.group())
        except json.JSONDecodeError:
            pass
    return None


# ============================================================
# A2: Design Test Plan
# ============================================================

async def a2_design_test(session, question: str, args) -> dict:
    """A2 designs a test plan based on the question (before seeing code)."""
    mol_list = "\n".join(f"  {i+1}. {s}" for i, s in enumerate(SINGLE_MOLECULES[:10]))

    messages = [
        {"role": "system", "content": A2_SYSTEM},
        {"role": "user", "content": A2_DESIGN_USER.format(
            question=question, molecules=mol_list,
        )},
    ]

    result = await call_llm(
        session, args.base_url, args.api_key,
        args.model_tester or args.model,
        messages, temperature=0.0, max_tokens=1024, timeout=args.timeout,
    )

    if "error" in result:
        return {"error": result["error"]}

    plan = parse_json_response(result["content"])
    if plan is None:
        # Fallback: use default test plan
        plan = {
            "test_smiles": [SINGLE_MOLECULES[0], SINGLE_MOLECULES[1]],
            "expected_type": "unknown",
            "value_description": "chemically reasonable output",
            "reasonable_range": "N/A",
            "failure_modes": ["runtime error", "wrong output type"],
        }

    # Validate test_smiles — must be from SINGLE_MOLECULES or valid SMILES
    if "test_smiles" not in plan or not plan["test_smiles"]:
        plan["test_smiles"] = [SINGLE_MOLECULES[0], SINGLE_MOLECULES[1]]

    return plan


# ============================================================
# A1: Generate Code
# ============================================================

async def a1_generate(session, question: str, args) -> dict:
    """A1 generates code for the question."""
    messages = [
        {"role": "system", "content": A1_SYSTEM},
        {"role": "user", "content": A1_GENERATE_USER.format(question=question)},
    ]

    result = await call_llm(
        session, args.base_url, args.api_key,
        args.model_coder or args.model,
        messages, temperature=args.temperature,
        max_tokens=args.max_tokens, timeout=args.timeout,
    )

    if "error" in result:
        return {"error": result["error"]}

    code = extract_code(result["content"])
    return {"code": code}


# ============================================================
# A1: Fix Code Based on Feedback
# ============================================================

async def a1_fix(session, question: str, code: str,
                 exec_results: str, feedback: str, args) -> dict:
    """A1 fixes code based on A2's feedback."""
    messages = [
        {"role": "system", "content": A1_SYSTEM},
        {"role": "user", "content": A1_FIX_USER.format(
            question=question, code=code,
            exec_results=exec_results, feedback=feedback,
        )},
    ]

    result = await call_llm(
        session, args.base_url, args.api_key,
        args.model_coder or args.model,
        messages, temperature=args.temperature,
        max_tokens=args.max_tokens, timeout=args.timeout,
    )

    if "error" in result:
        return {"error": result["error"]}

    fixed = extract_code(result["content"])
    return {"code": fixed}


# ============================================================
# A2: Verify Code
# ============================================================

async def a2_verify(session, question: str, test_plan: dict,
                    exec_results: str, args) -> dict:
    """A2 judges whether the execution results are chemically correct."""
    messages = [
        {"role": "system", "content": A2_SYSTEM},
        {"role": "user", "content": A2_VERIFY_USER.format(
            question=question,
            test_plan=json.dumps(test_plan, indent=2),
            exec_results=exec_results,
            expected_type=test_plan.get("expected_type", "unknown"),
            value_description=test_plan.get("value_description", ""),
        )},
    ]

    result = await call_llm(
        session, args.base_url, args.api_key,
        args.model_tester or args.model,
        messages, temperature=0.0, max_tokens=1024, timeout=args.timeout,
    )

    if "error" in result:
        return {"verdict": "ERROR", "reason": result["error"]}

    verdict = parse_json_response(result["content"])
    if verdict is None:
        # If we can't parse, treat as PASS (conservative — don't block)
        return {"verdict": "PASS", "reason": "unparseable response, defaulting to PASS"}

    return verdict


# ============================================================
# Execute Code with Multiple Test Molecules
# ============================================================

def run_tests(code: str, test_smiles_list: list, exec_timeout: int = 30) -> str:
    """Execute code with each test molecule and collect results."""
    results = []
    for smi in test_smiles_list:
        r = execute_code(code, smi, timeout=exec_timeout)
        results.append(f"SMILES: {smi}\n  Result: {r}")
    return "\n".join(results)


# ============================================================
# Single Question: Dual Agent Loop
# ============================================================

async def run_single_dual(
    semaphore: asyncio.Semaphore,
    session,
    level: int,
    idx: int,
    question: str,
    output_dir: Path,
    args,
) -> dict:
    """Dual agent loop for a single question."""
    output_path = output_dir / f"temp{idx}.py"

    if args.resume and output_path.exists():
        return {"level": level, "idx": idx, "status": "skipped",
                "rounds": 0, "api_calls": 0}

    async with semaphore:
        api_calls = 0

        # Step 1: A2 designs test plan
        test_plan = await a2_design_test(session, question, args)
        api_calls += 1
        if "error" in test_plan:
            return {"level": level, "idx": idx, "status": "api_error",
                    "error": test_plan["error"], "rounds": 0, "api_calls": api_calls}

        test_smiles = test_plan.get("test_smiles", [SINGLE_MOLECULES[0]])
        # Limit to 3 test molecules max
        test_smiles = test_smiles[:3]

        # Step 2: A1 generates code
        gen_result = await a1_generate(session, question, args)
        api_calls += 1
        if "error" in gen_result:
            return {"level": level, "idx": idx, "status": "api_error",
                    "error": gen_result["error"], "rounds": 0, "api_calls": api_calls}

        code = gen_result["code"]
        if not code or "level_function" not in code:
            return {"level": level, "idx": idx, "status": "failed",
                    "rounds": 0, "api_calls": api_calls}

        # Step 3-6: Test-Fix loop
        final_code = code
        passed = False

        for round_num in range(1, args.max_rounds + 1):
            # Execute code with test molecules
            exec_results = run_tests(final_code, test_smiles, args.exec_timeout)

            # A2 verifies
            verdict = await a2_verify(session, question, test_plan, exec_results, args)
            api_calls += 1

            if verdict.get("verdict", "").upper() == "PASS":
                passed = True
                break

            # Last round — don't try to fix, just save what we have
            if round_num >= args.max_rounds:
                break

            # A1 fixes based on feedback
            feedback = verdict.get("reason", "Test failed, please review and fix.")
            fix_result = await a1_fix(
                session, question, final_code, exec_results, feedback, args
            )
            api_calls += 1

            if "error" in fix_result:
                break  # API error during fix, save current code

            fixed_code = fix_result["code"]
            if fixed_code and "level_function" in fixed_code:
                final_code = fixed_code
            else:
                break  # Fix produced invalid code, save current

        # Save final code
        if final_code and "level_function" in final_code:
            output_dir.mkdir(parents=True, exist_ok=True)
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(final_code)
            return {
                "level": level, "idx": idx,
                "status": "success",
                "passed_test": passed,
                "rounds": round_num if 'round_num' in dir() else 0,
                "api_calls": api_calls,
            }
        else:
            return {
                "level": level, "idx": idx,
                "status": "failed",
                "rounds": round_num if 'round_num' in dir() else 0,
                "api_calls": api_calls,
            }


# ============================================================
# Main Runner
# ============================================================

async def run_dual_agent(args):
    import aiohttp

    model_tag = args.model.replace("/", "_").replace(":", "_")
    run_tag = f"{model_tag}_{args.lang}_dual"
    if args.tag:
        run_tag = args.tag

    base_output = PROJECT_ROOT / "predictions" / run_tag

    print("=" * 60)
    print("  MolViBench — Dual Agent (Coder-Tester) Inference")
    print("=" * 60)
    print(f"  Model:       {args.model}")
    if args.model_coder:
        print(f"  Coder model: {args.model_coder}")
    if args.model_tester:
        print(f"  Tester model:{args.model_tester}")
    print(f"  Base URL:    {args.base_url}")
    print(f"  Language:    {args.lang}")
    print(f"  Max rounds:  {args.max_rounds}")
    print(f"  Output:      {base_output}")
    print(f"  Concurrency: {args.concurrency}")
    print(f"  Resume:      {args.resume}")
    print("=" * 60)

    tasks_info = []
    for level in args.levels:
        csv_path = PROJECT_ROOT / "data" / args.lang / f"level{level}.csv"
        if not csv_path.exists():
            print(f"  [WARN] CSV not found: {csv_path}")
            continue
        questions = load_questions(str(csv_path))
        level_dir = base_output / f"level{level}"
        for i, q in enumerate(questions, 1):
            tasks_info.append((level, i, q, level_dir))

    total = len(tasks_info)
    print(f"\n  Total tasks: {total}")

    if args.resume:
        already = sum(1 for _, i, _, ld in tasks_info if (ld / f"temp{i}.py").exists())
        print(f"  Already done: {already}")
        print(f"  Remaining: {total - already}")

    semaphore = asyncio.Semaphore(args.concurrency)
    start_time = time.time()

    stats = {
        "success": 0, "failed": 0, "api_error": 0, "skipped": 0,
        "passed_test": 0, "total_rounds": 0, "total_api_calls": 0,
    }
    done_count = 0

    async def wrapped(level, idx, question, output_dir):
        nonlocal done_count
        r = await run_single_dual(semaphore, session, level, idx, question, output_dir, args)
        done_count += 1

        if r["status"] == "success":
            stats["success"] += 1
            if r.get("passed_test"):
                stats["passed_test"] += 1
            stats["total_rounds"] += r.get("rounds", 0)
            stats["total_api_calls"] += r.get("api_calls", 0)
        elif r["status"] == "failed":
            stats["failed"] += 1
            stats["total_rounds"] += r.get("rounds", 0)
            stats["total_api_calls"] += r.get("api_calls", 0)
        elif r["status"] == "api_error":
            stats["api_error"] += 1
        elif r["status"] == "skipped":
            stats["skipped"] += 1

        if done_count % 5 == 0 or done_count == total:
            elapsed = time.time() - start_time
            rate = elapsed / max(done_count - stats["skipped"], 1)
            remaining = total - done_count
            eta = int(rate * remaining)
            ok = stats["success"]
            fail = stats["failed"]
            err = stats["api_error"]
            skip = stats["skipped"]
            pt = stats["passed_test"]
            print(f"  [{done_count:>4}/{total}] ok={ok}(pass={pt}) fail={fail} "
                  f"err={err} skip={skip} | ETA {eta}s")
        return r

    connector = aiohttp.TCPConnector(limit=args.concurrency + 2, ssl=False)
    async with aiohttp.ClientSession(connector=connector) as session:
        coros = [wrapped(lv, idx, q, od) for lv, idx, q, od in tasks_info]
        results = await asyncio.gather(*coros)

    elapsed = time.time() - start_time

    print(f"\n{'=' * 60}")
    print("  DUAL AGENT INFERENCE COMPLETE")
    print("=" * 60)
    print(f"  Total:       {total}")
    print(f"  Success:     {stats['success']}")
    print(f"  Passed test: {stats['passed_test']}")
    print(f"  Failed:      {stats['failed']}")
    print(f"  API Errors:  {stats['api_error']}")
    print(f"  Skipped:     {stats['skipped']}")
    print(f"  Time:        {elapsed:.1f}s")

    n_active = stats["success"] + stats["failed"]
    if n_active > 0:
        avg_rounds = stats["total_rounds"] / n_active
        avg_api = stats["total_api_calls"] / n_active
        print(f"\n  Avg rounds/question:     {avg_rounds:.1f}")
        print(f"  Avg API calls/question:  {avg_api:.1f}")
        print(f"  Test pass rate:          {stats['passed_test']}/{n_active} "
              f"({stats['passed_test']/n_active:.1%})")

    # API errors
    api_errors = [r for r in results if r["status"] == "api_error"]
    if api_errors:
        print(f"\n  First 5 API errors:")
        for r in api_errors[:5]:
            print(f"    L{r['level']}/temp{r['idx']}: {str(r.get('error', ''))[:80]}")

    # Save metadata
    meta = {
        "model": args.model,
        "model_coder": args.model_coder,
        "model_tester": args.model_tester,
        "base_url": args.base_url,
        "lang": args.lang,
        "mode": "dual_agent",
        "max_rounds": args.max_rounds,
        "temperature": args.temperature,
        "max_tokens": args.max_tokens,
        "levels": args.levels,
        "concurrency": args.concurrency,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "elapsed_s": round(elapsed, 1),
        "stats": stats,
    }
    meta_path = base_output / "run_meta.json"
    base_output.mkdir(parents=True, exist_ok=True)
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
    print(f"  Metadata:    {meta_path}")
    print(f"  Output:      {base_output}")


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="MolViBench Dual Agent (Coder-Tester) Inference"
    )
    parser.add_argument("--base-url", type=str, required=True)
    parser.add_argument("--api-key", type=str, required=True)
    parser.add_argument("--model", type=str, required=True,
                        help="Default model for both agents")
    parser.add_argument("--model-coder", type=str, default=None,
                        help="Override model for A1 (Coder)")
    parser.add_argument("--model-tester", type=str, default=None,
                        help="Override model for A2 (Tester)")
    parser.add_argument("--lang", type=str, default="en", choices=["en", "cn"])
    parser.add_argument("--levels", type=int, nargs="+", default=[1, 2, 3, 4, 5])
    parser.add_argument("--max-rounds", type=int, default=3,
                        help="Max A1-A2 interaction rounds per question")
    parser.add_argument("--temperature", type=float, default=0.0)
    parser.add_argument("--max-tokens", type=int, default=4096)
    parser.add_argument("--timeout", type=int, default=120,
                        help="API call timeout (seconds)")
    parser.add_argument("--exec-timeout", type=int, default=30,
                        help="Code execution timeout")
    parser.add_argument("--concurrency", type=int, default=3)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--tag", type=str, default=None)

    args = parser.parse_args()
    asyncio.run(run_dual_agent(args))


if __name__ == "__main__":
    main()
