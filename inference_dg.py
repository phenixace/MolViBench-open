#!/usr/bin/env python
"""
MolViBench - Direct Generation (DG)
=====================================

Single-pass code generation: LLM receives task description and produces code directly.
No execution feedback or iterative repair.

Usage:
    python inference_dg.py --base-url https://api.openai.com/v1 --api-key sk-xxx --model gpt-4o
    python inference_dg.py --base-url http://localhost:8000/v1 --api-key EMPTY --model my-model
    python inference_dg.py ... --levels 1 2 3 --lang en --prompt-strategy signature --resume
"""

import argparse
import asyncio
import csv
import json
import os
import sys
import time
import re
from pathlib import Path
from typing import Optional

# Fix Windows console encoding and asyncio event loop
if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())

PROJECT_ROOT = Path(__file__).resolve().parent

# ============================================================
# Prompt construction
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


def build_user_prompt(question: str, lang: str, strategy: str) -> str:
    """Build user prompt from question text and strategy."""
    if strategy == "basic":
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
    else:
        raise ValueError(f"Unknown prompt strategy: {strategy}")


def get_system_prompt(lang: str) -> str:
    return SYSTEM_PROMPT_CN if lang == "cn" else SYSTEM_PROMPT_EN


# ============================================================
# Code extraction
# ============================================================

def extract_code(response: str) -> str:
    """Extract Python code from LLM response (handles markdown blocks)."""
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
# Question loader
# ============================================================

def load_questions(csv_path: str) -> list[str]:
    """Load questions from CSV. Returns list of question strings."""
    questions = []
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            questions.append(row.get("question", "").strip())
    return questions


# ============================================================
# Async API caller
# ============================================================

async def call_api(
    session,
    base_url: str,
    api_key: str,
    model: str,
    system_prompt: str,
    user_prompt: str,
    temperature: float = 0.0,
    max_tokens: int = 4096,
    timeout: int = 120,
) -> dict:
    """Call OpenAI-compatible chat completion API."""
    import aiohttp
    import ssl as _ssl

    url = f"{base_url.rstrip('/')}/chat/completions"
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {api_key}",
    }
    payload = {
        "model": model,
        "messages": [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ],
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
# Task: one question
# ============================================================

def _is_retriable_error(error: str) -> bool:
    """True if error is timeout or HTTP 524 (retry may help)."""
    if not error:
        return False
    err_lower = error.lower()
    return "timeout" in err_lower or "524" in err_lower or "request timeout" in err_lower


async def process_question(
    semaphore: asyncio.Semaphore,
    session,
    level: int,
    idx: int,
    question: str,
    output_dir: Path,
    args,
) -> dict:
    """Process a single question: build prompt → call API (with retries on timeout) → save code."""
    output_path = output_dir / f"temp{idx}.py"

    if args.resume and output_path.exists():
        return {
            "level": level,
            "idx": idx,
            "status": "skipped",
            "reason": "already exists",
        }

    system_prompt = get_system_prompt(args.lang)
    user_prompt = build_user_prompt(question, args.lang, args.prompt_strategy)
    max_retries = getattr(args, "max_retries", 1)
    last_result = None

    for attempt in range(max_retries):
        async with semaphore:
            result = await call_api(
                session=session,
                base_url=args.base_url,
                api_key=args.api_key,
                model=args.model,
                system_prompt=system_prompt,
                user_prompt=user_prompt,
                temperature=args.temperature,
                max_tokens=args.max_tokens,
                timeout=args.timeout,
            )
        last_result = result
        if result["success"]:
            break
        if not _is_retriable_error(result.get("error", "")) or attempt == max_retries - 1:
            break
        # retriable and still have attempts left
        await asyncio.sleep(2.0 * (attempt + 1))  # brief backoff before retry

    if not last_result or not last_result["success"]:
        return {
            "level": level,
            "idx": idx,
            "status": "error",
            "error": last_result["error"] if last_result else "No response",
        }

    code = extract_code(last_result["content"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(code, encoding="utf-8")

    return {
        "level": level,
        "idx": idx,
        "status": "success",
        "usage": last_result.get("usage", {}),
        "code_length": len(code),
    }


# ============================================================
# Main inference loop
# ============================================================

async def run_inference(args):
    """Run inference across all specified levels."""
    import aiohttp

    sys.stdout.flush()

    # Determine output directory
    model_tag = args.model.replace("/", "_").replace(":", "_")
    run_tag = f"{model_tag}_{args.lang}_{args.prompt_strategy}"
    if args.tag:
        run_tag = args.tag

    base_output = PROJECT_ROOT / "predictions" / run_tag

    print("=" * 60)
    print("  MolViBench — Direct LLM Inference")
    print("=" * 60)
    print(f"  Model:     {args.model}")
    print(f"  Base URL:  {args.base_url}")
    print(f"  Language:  {args.lang}")
    print(f"  Strategy:  {args.prompt_strategy}")
    print(f"  Output:    {base_output}")
    print(f"  Concurrency: {args.concurrency}")
    print(f"  Resume:    {args.resume}")
    print(f"  Max retries (on timeout/524): {getattr(args, 'max_retries', 1)}")
    print(f"  Temperature: {args.temperature}")
    print("=" * 60)

    # Collect all tasks
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

    # Count already done (for resume)
    already_done = 0
    if args.resume:
        for level, idx, q, output_dir in tasks_info:
            if (output_dir / f"temp{idx}.py").exists():
                already_done += 1
        print(f"  Already completed: {already_done}")
        print(f"  Remaining: {total - already_done}")

    semaphore = asyncio.Semaphore(args.concurrency)
    start_time = time.time()

    stats = {"success": 0, "error": 0, "skipped": 0, "total_tokens": 0}
    errors = []

    connector = aiohttp.TCPConnector(limit=args.concurrency * 2)
    async with aiohttp.ClientSession(connector=connector) as session:
        coros = [
            process_question(semaphore, session, level, idx, q, output_dir, args)
            for level, idx, q, output_dir in tasks_info
        ]

        done_count = 0
        for coro in asyncio.as_completed(coros):
            result = await coro
            done_count += 1

            status = result["status"]
            stats[status] = stats.get(status, 0) + 1

            if status == "success":
                usage = result.get("usage", {})
                stats["total_tokens"] += usage.get("total_tokens", 0)

            if status == "error":
                errors.append(result)

            # Progress
            if done_count % 10 == 0 or done_count == total:
                elapsed = time.time() - start_time
                rate = done_count / elapsed if elapsed > 0 else 0
                eta = (total - done_count) / rate if rate > 0 else 0
                print(
                    f"  [{done_count:4d}/{total}] "
                    f"ok={stats['success']} err={stats['error']} skip={stats['skipped']} "
                    f"| {rate:.1f} q/s | ETA {eta:.0f}s",
                    flush=True,
                )

    elapsed = time.time() - start_time

    # Summary
    print(f"\n{'='*60}")
    print("  INFERENCE COMPLETE")
    print(f"{'='*60}")
    print(f"  Total:    {total}")
    print(f"  Success:  {stats['success']}")
    print(f"  Errors:   {stats['error']}")
    print(f"  Skipped:  {stats['skipped']}")
    print(f"  Tokens:   {stats['total_tokens']:,}")
    print(f"  Time:     {elapsed:.1f}s")
    print(f"  Output:   {base_output}")

    if errors:
        print(f"\n  First 5 errors:")
        for e in errors[:5]:
            print(f"    L{e['level']}/temp{e['idx']}: {e['error'][:100]}")

    # Save run metadata
    meta = {
        "model": args.model,
        "base_url": args.base_url,
        "lang": args.lang,
        "prompt_strategy": args.prompt_strategy,
        "temperature": args.temperature,
        "max_tokens": args.max_tokens,
        "levels": args.levels,
        "concurrency": args.concurrency,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "elapsed_s": round(elapsed, 1),
        "stats": stats,
    }
    meta_path = base_output / "run_meta.json"
    meta_path.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"  Metadata: {meta_path}")

    return stats


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="MolViBench — Direct LLM Inference",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # API config
    parser.add_argument("--base-url", required=True, help="OpenAI-compatible API base URL")
    parser.add_argument("--api-key", required=True, help="API key")
    parser.add_argument("--model", required=True, help="Model name (e.g. gpt-4o, deepseek-chat)")

    # Task config
    parser.add_argument(
        "--levels", type=int, nargs="+", default=[1, 2, 3, 4, 5],
        help="Levels to evaluate (default: 1 2 3 4 5)",
    )
    parser.add_argument("--lang", default="en", choices=["cn", "en"], help="Prompt language")
    parser.add_argument(
        "--prompt-strategy", default="basic", choices=["basic"],
        help="Prompt strategy",
    )

    # Generation config
    parser.add_argument("--temperature", type=float, default=0.0, help="Sampling temperature")
    parser.add_argument("--max-tokens", type=int, default=4096, help="Max tokens per response")
    parser.add_argument("--timeout", type=int, default=120, help="API request timeout (seconds)")

    # Execution config
    parser.add_argument("--concurrency", type=int, default=10, help="Max concurrent API requests")
    parser.add_argument("--resume", action="store_true", help="Skip already-generated predictions")
    parser.add_argument("--tag", type=str, default=None, help="Custom run tag (output folder name)")
    parser.add_argument(
        "--max-retries", type=int, default=1,
        help="Max retries per question on timeout/524 (default 1 = no retry)",
    )

    args = parser.parse_args()
    asyncio.run(run_inference(args))


if __name__ == "__main__":
    main()
