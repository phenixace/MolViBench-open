#!/usr/bin/env python
"""
MolViBench - Incremental Repair (IR)
===============================

Iterative Incremental Repair: LLM generates code, code is executed locally,
if execution fails the error is fed back to the LLM for correction.
Repeats up to N rounds.

Usage:
    python inference_ir.py --base-url https://api.openai.com/v1 --api-key sk-xxx --model gpt-4o
    python inference_ir.py ... --max-rounds 3 --lang cn --resume
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
import traceback
from pathlib import Path
from typing import Optional

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from evaluate.molecules import SINGLE_MOLECULES, REPRESENTATIVE_MOLECULES


# ============================================================
# Prompts
# ============================================================

SYSTEM_PROMPT_EN = """You are an expert cheminformatics programmer.
Write Python code using RDKit to solve molecular computing tasks.
Your code must define a function called `level_function` that takes the specified input and returns the result.
Only use these libraries: rdkit, numpy, pandas, scikit-learn, matplotlib, selfies.
Include error handling with try/except blocks.
Do NOT include any import statements for libraries other than the ones listed above.
Return only the Python code, no explanations."""

SYSTEM_PROMPT_CN = """你是一名化学信息学编程专家。
请使用 RDKit 编写 Python 代码解决分子计算任务。
你的代码必须定义一个名为 `level_function` 的函数，接受指定输入并返回结果。
仅使用以下库：rdkit, numpy, pandas, scikit-learn, matplotlib, selfies。
包含 try/except 错误处理。
不要导入上述列表之外的库。
只返回 Python 代码，不要解释。"""

REPAIR_PROMPT_EN = """The code you provided failed during execution. Here is the error:

```
{error}
```

The test was run as:
```python
result = level_function({test_input})
```

Please fix the code and provide the corrected complete Python code.
The function must still be named `level_function`.
Return only the Python code, no explanations."""

REPAIR_PROMPT_CN = """你提供的代码在执行时失败了。以下是错误信息：

```
{error}
```

测试调用方式：
```python
result = level_function({test_input})
```

请修复代码并给出完整的修正后的 Python 代码。
函数名必须仍为 `level_function`。
只返回 Python 代码，不要解释。"""


def get_system_prompt(lang: str) -> str:
    return SYSTEM_PROMPT_CN if lang == "cn" else SYSTEM_PROMPT_EN


def build_initial_prompt(question: str, lang: str) -> str:
    if lang == "cn":
        return (
            f"请编写一个 Python 函数 `level_function` 来完成以下任务：\n\n"
            f"任务：{question}\n\n"
            f"要求：\n"
            f"1. 函数名必须为 `level_function`\n"
            f"2. 使用 RDKit 库\n"
            f"3. 包含 try/except 错误处理\n"
            f"4. 输入无效时返回 None\n\n"
            f"请直接给出完整的 Python 代码。"
        )
    else:
        return (
            f"Write a Python function `level_function` to accomplish the following task:\n\n"
            f"Task: {question}\n\n"
            f"Requirements:\n"
            f"1. Function must be named `level_function`\n"
            f"2. Use the RDKit library\n"
            f"3. Include try/except error handling\n"
            f"4. Return None for invalid input\n\n"
            f"Provide the complete Python code."
        )


def build_repair_prompt(error: str, test_input: str, lang: str) -> str:
    template = REPAIR_PROMPT_CN if lang == "cn" else REPAIR_PROMPT_EN
    return template.format(error=error, test_input=test_input)


# ============================================================
# Code extraction
# ============================================================

def extract_code(response: str) -> str:
    """Extract Python code from LLM response."""
    patterns = [
        r"```python\s*\n(.*?)```",
        r"```py\s*\n(.*?)```",
        r"```\s*\n(.*?)```",
    ]
    for pattern in patterns:
        matches = re.findall(pattern, response, re.DOTALL)
        if matches:
            return max(matches, key=len).strip()

    lines = response.strip().split("\n")
    if any("def " in l for l in lines) and any(
        "import" in l or "from " in l for l in lines
    ):
        return response.strip()
    return response.strip()


# ============================================================
# Local code execution (for Incremental Repair feedback)
# ============================================================

def execute_code_locally(code: str, test_smiles: str, timeout: int = 30) -> dict:
    """
    Execute generated code with a test SMILES in a subprocess.
    Writes the user code to a temp file and a separate runner script
    to avoid string-escaping issues with exec().
    """
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
    result = {{"success": False, "error": f"Code load error: {{type(e).__name__}}: {{e}}", "output": None}}
    print("__AGENT_RESULT__" + json.dumps(result, default=str))
    sys.exit(0)

if not hasattr(mod, "level_function"):
    result = {{"success": False, "error": "Function 'level_function' not defined in the code.", "output": None}}
    print("__AGENT_RESULT__" + json.dumps(result, default=str))
    sys.exit(0)

try:
    output = mod.level_function({repr(test_smiles)})
    result = {{"success": True, "error": None, "output": repr(output), "output_type": type(output).__name__}}
    print("__AGENT_RESULT__" + json.dumps(result, default=str))
except Exception as e:
    tb = traceback.format_exc()
    result = {{"success": False, "error": f"{{type(e).__name__}}: {{e}}\\n{{tb}}", "output": None}}
    print("__AGENT_RESULT__" + json.dumps(result, default=str))
'''

    with open(runner_path, "w", encoding="utf-8") as f:
        f.write(runner)

    try:
        proc = subprocess.run(
            [sys.executable, runner_path],
            capture_output=True,
            text=True,
            timeout=timeout,
            env={**os.environ, "PYTHONIOENCODING": "utf-8"},
        )

        stdout = proc.stdout
        stderr = proc.stderr

        for line in stdout.split("\n"):
            if line.startswith("__AGENT_RESULT__"):
                data = json.loads(line[len("__AGENT_RESULT__"):])
                return data

        return {
            "success": False,
            "error": f"No result returned.\nstdout: {stdout[:500]}\nstderr: {stderr[:500]}",
            "output": None,
        }

    except subprocess.TimeoutExpired:
        return {
            "success": False,
            "error": f"Execution timeout ({timeout}s)",
            "output": None,
        }
    except Exception as e:
        return {"success": False, "error": str(e), "output": None}
    finally:
        import shutil
        try:
            shutil.rmtree(code_dir, ignore_errors=True)
        except OSError:
            pass


# ============================================================
# API caller
# ============================================================

async def call_api(
    session,
    base_url: str,
    api_key: str,
    model: str,
    messages: list[dict],
    temperature: float = 0.0,
    max_tokens: int = 4096,
    timeout: int = 120,
) -> dict:
    """Call OpenAI-compatible chat completion API with full message history."""
    import aiohttp
    import ssl as _ssl

    url = f"{base_url.rstrip('/')}/chat/completions"
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {api_key}",
    }
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
            url,
            headers=headers,
            json=payload,
            timeout=aiohttp.ClientTimeout(total=timeout),
            ssl=ssl_ctx,
        ) as resp:
            if resp.status != 200:
                text = await resp.text()
                return {
                    "success": False,
                    "error": f"HTTP {resp.status}: {text[:500]}",
                    "content": "",
                }
            data = await resp.json()
            content = data["choices"][0]["message"]["content"]
            usage = data.get("usage", {})
            return {
                "success": True,
                "content": content,
                "usage": usage,
                "error": None,
            }
    except asyncio.TimeoutError:
        return {"success": False, "error": "Request timeout", "content": ""}
    except Exception as e:
        return {"success": False, "error": str(e), "content": ""}


# ============================================================
# Multi-molecule execution test
# ============================================================

def execute_code_multi_molecule(code: str, molecules: list, timeout: int = 30) -> dict:
    """Test code against multiple representative molecules.
    Returns success=True if ANY molecule executes without error.
    """
    last_error = "No molecules tested"
    for smi in molecules:
        result = execute_code_locally(code, smi, timeout=timeout)
        if result["success"]:
            return result
        last_error = result["error"]
    return {"success": False, "error": last_error, "output": None}


# ============================================================
# Agent: one question with Incremental Repair loop
# ============================================================

async def agent_solve_question(
    semaphore: asyncio.Semaphore,
    session,
    level: int,
    idx: int,
    question: str,
    output_dir: Path,
    args,
) -> dict:
    """
    Agent loop for a single question:
    1. Generate initial code
    2. Execute locally
    3. If error → feed back and retry
    4. Save the best/final code
    """
    output_path = output_dir / f"temp{idx}.py"

    if args.resume and output_path.exists():
        return {
            "level": level, "idx": idx, "status": "skipped",
            "rounds": 0, "reason": "already exists",
        }

    test_smiles = REPRESENTATIVE_MOLECULES[0]  # for repair prompt context
    test_molecules = REPRESENTATIVE_MOLECULES   # for execution testing
    system_prompt = get_system_prompt(args.lang)
    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": build_initial_prompt(question, args.lang)},
    ]

    best_code = None
    final_round = 0

    async with semaphore:
        for round_num in range(1, args.max_rounds + 1):
            # Call LLM
            result = await call_api(
                session=session,
                base_url=args.base_url,
                api_key=args.api_key,
                model=args.model,
                messages=messages,
                temperature=args.temperature,
                max_tokens=args.max_tokens,
                timeout=args.timeout,
            )

            if not result["success"]:
                return {
                    "level": level, "idx": idx, "status": "api_error",
                    "rounds": round_num, "error": result["error"],
                }

            code = extract_code(result["content"])
            best_code = code
            final_round = round_num

            # Execute locally with representative molecules
            exec_result = execute_code_multi_molecule(
                code, test_molecules, timeout=args.exec_timeout
            )

            if exec_result["success"]:
                break

            # Execution failed — prepare repair prompt
            error_msg = exec_result["error"]
            if len(error_msg) > 1500:
                error_msg = error_msg[:1500] + "\n... (truncated)"

            repair_prompt = build_repair_prompt(
                error=error_msg,
                test_input=repr(test_smiles),
                lang=args.lang,
            )

            messages.append({"role": "assistant", "content": result["content"]})
            messages.append({"role": "user", "content": repair_prompt})

    # Save final code
    if best_code:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(best_code, encoding="utf-8")

    status = "success" if exec_result["success"] else "failed_after_repair"
    return {
        "level": level,
        "idx": idx,
        "status": status,
        "rounds": final_round,
        "output_type": exec_result.get("output_type", "unknown"),
    }


# ============================================================
# Question loader
# ============================================================

def load_questions(csv_path: str) -> list[str]:
    questions = []
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            questions.append(row.get("question", "").strip())
    return questions


# ============================================================
# Main
# ============================================================

async def run_agent(args):
    import aiohttp

    model_tag = args.model.replace("/", "_").replace(":", "_")
    run_tag = f"{model_tag}_{args.lang}_agent_r{args.max_rounds}"
    if args.tag:
        run_tag = args.tag

    base_output = PROJECT_ROOT / "predictions" / run_tag

    print("=" * 60)
    print("  MolViBench — Agent Mode Inference (Incremental Repair)")
    print("=" * 60)
    print(f"  Model:       {args.model}")
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
        output_dir = base_output / f"level{level}"
        output_dir.mkdir(parents=True, exist_ok=True)
        for idx, q in enumerate(questions, 1):
            tasks_info.append((level, idx, q, output_dir))

    total = len(tasks_info)
    print(f"\n  Total tasks: {total}")

    if args.resume:
        already = sum(
            1 for l, i, q, d in tasks_info if (d / f"temp{i}.py").exists()
        )
        print(f"  Already completed: {already}")
        print(f"  Remaining: {total - already}")

    semaphore = asyncio.Semaphore(args.concurrency)
    start_time = time.time()

    stats = {
        "success": 0, "failed_after_repair": 0, "api_error": 0, "skipped": 0,
        "rounds_distribution": {},
    }
    errors = []

    connector = aiohttp.TCPConnector(limit=args.concurrency * 2)
    async with aiohttp.ClientSession(connector=connector) as session:
        coros = [
            agent_solve_question(
                semaphore, session, level, idx, q, output_dir, args
            )
            for level, idx, q, output_dir in tasks_info
        ]

        done_count = 0
        for coro in asyncio.as_completed(coros):
            result = await coro
            done_count += 1

            status = result["status"]
            stats[status] = stats.get(status, 0) + 1

            if status in ("success", "failed_after_repair"):
                r = result["rounds"]
                stats["rounds_distribution"][r] = (
                    stats["rounds_distribution"].get(r, 0) + 1
                )

            if status == "api_error":
                errors.append(result)

            # Progress
            if done_count % 5 == 0 or done_count == total:
                elapsed = time.time() - start_time
                rate = done_count / elapsed if elapsed > 0 else 0
                eta = (total - done_count) / rate if rate > 0 else 0
                rnd_info = " ".join(
                    f"R{k}={v}" for k, v in sorted(stats["rounds_distribution"].items())
                )
                print(
                    f"  [{done_count:4d}/{total}] "
                    f"ok={stats['success']} fail={stats['failed_after_repair']} "
                    f"err={stats['api_error']} skip={stats['skipped']} "
                    f"| {rnd_info} "
                    f"| ETA {eta:.0f}s",
                    flush=True,
                )

    elapsed = time.time() - start_time

    print(f"\n{'='*60}")
    print("  AGENT INFERENCE COMPLETE")
    print(f"{'='*60}")
    print(f"  Total:       {total}")
    print(f"  Success:     {stats['success']}")
    print(f"  Failed:      {stats['failed_after_repair']}")
    print(f"  API Errors:  {stats['api_error']}")
    print(f"  Skipped:     {stats['skipped']}")
    print(f"  Time:        {elapsed:.1f}s")
    print(f"\n  Rounds distribution (how many rounds needed):")
    for r in sorted(stats["rounds_distribution"]):
        cnt = stats["rounds_distribution"][r]
        print(f"    Round {r}: {cnt} questions")

    if errors:
        print(f"\n  First 5 API errors:")
        for e in errors[:5]:
            print(f"    L{e['level']}/temp{e['idx']}: {e.get('error', '')[:100]}")

    # Save metadata
    meta = {
        "model": args.model,
        "base_url": args.base_url,
        "lang": args.lang,
        "mode": "agent",
        "max_rounds": args.max_rounds,
        "temperature": args.temperature,
        "max_tokens": args.max_tokens,
        "levels": args.levels,
        "concurrency": args.concurrency,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "elapsed_s": round(elapsed, 1),
        "stats": {k: v for k, v in stats.items()},
    }
    meta_path = base_output / "run_meta.json"
    meta_path.write_text(
        json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"  Metadata:    {meta_path}")
    print(f"  Output:      {base_output}")


def main():
    parser = argparse.ArgumentParser(
        description="MolViBench — Agent Mode Inference (Incremental Repair)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("--base-url", required=True, help="OpenAI-compatible API base URL")
    parser.add_argument("--api-key", required=True, help="API key")
    parser.add_argument("--model", required=True, help="Model name")

    parser.add_argument(
        "--levels", type=int, nargs="+", default=[1, 2, 3, 4, 5],
        help="Levels to evaluate (default: 1 2 3 4 5)",
    )
    parser.add_argument("--lang", default="en", choices=["cn", "en"])
    parser.add_argument("--max-rounds", type=int, default=3, help="Max Incremental Repair rounds per question")

    parser.add_argument("--temperature", type=float, default=0.0)
    parser.add_argument("--max-tokens", type=int, default=4096)
    parser.add_argument("--timeout", type=int, default=120, help="API request timeout")
    parser.add_argument("--exec-timeout", type=int, default=30, help="Local code execution timeout")

    parser.add_argument("--concurrency", type=int, default=5, help="Max concurrent questions")
    parser.add_argument("--resume", action="store_true", help="Skip already-generated predictions")
    parser.add_argument("--tag", type=str, default=None, help="Custom run tag (output folder name)")

    args = parser.parse_args()
    asyncio.run(run_agent(args))


if __name__ == "__main__":
    main()
