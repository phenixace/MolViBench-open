#!/usr/bin/env python
"""
RDKitBench — Molecular Vibe Coding Benchmark
=============================================

Main entry point for evaluation.

Usage:
    # Evaluate all levels (LLM predictions vs solutions)
    python test.py --all

    # Evaluate a specific level
    python test.py --level 1

    # Evaluate a single question
    python test.py --level 1 --question 5

    # Generate ground truth cache
    python test.py --generate-ground-truth

    # Validate all solutions (self-test)
    python test.py --validate-solutions

    # Change language (cn/en)
    python test.py --all --lang en

    # Use direct execution (faster, less safe)
    python test.py --all --no-safe
"""

import argparse
import os
import sys
import json
import time

# Fix Windows console encoding
if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    sys.stderr.reconfigure(encoding='utf-8', errors='replace')

# Add project root to path
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PROJECT_ROOT)

from evaluate.evaluator import (
    EvalConfig, evaluate_all, evaluate_level, evaluate_single,
    save_report, load_questions, infer_test_inputs,
)
from evaluate.executor import (
    execute_function_direct, check_syntax, check_function_exists,
)
from evaluate.comparators import compare_value


def cmd_evaluate_all(args):
    """Evaluate all levels."""
    config = EvalConfig(PROJECT_ROOT)
    config.lang = args.lang
    config.safe_mode = not args.no_safe
    config.timeout = args.timeout
    config.n_test_molecules = args.n_test

    report = evaluate_all(config)

    # Save report
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    report_path = os.path.join(PROJECT_ROOT, "results", f"report_{args.lang}_{timestamp}.json")
    save_report(report, report_path)


def cmd_evaluate_level(args):
    """Evaluate a specific level."""
    config = EvalConfig(PROJECT_ROOT)
    config.lang = args.lang
    config.safe_mode = not args.no_safe
    config.timeout = args.timeout
    config.n_test_molecules = args.n_test
    config.levels = [args.level]

    if args.question:
        # Single question evaluation
        sol_path = config.solution_path(args.level, args.question)
        pred_path = config.prediction_path(args.level, args.question)
        result = evaluate_single(sol_path, pred_path, config)
        result["question_idx"] = args.question

        questions = load_questions(config.csv_path(args.level))
        if args.question <= len(questions):
            result["question"] = questions[args.question - 1]

        print(json.dumps(result, indent=2, ensure_ascii=False, default=str))
    else:
        summary = evaluate_level(args.level, config)
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        report_path = os.path.join(PROJECT_ROOT, "results",
                                   f"report_L{args.level}_{args.lang}_{timestamp}.json")
        save_report({"level_details": {args.level: summary}}, report_path)


def cmd_validate_solutions(args):
    """Self-test: validate all solution files can execute."""
    config = EvalConfig(PROJECT_ROOT)
    config.lang = args.lang

    print("=" * 60)
    print("  Validating all solution files")
    print("=" * 60)

    stats = {"total": 0, "syntax_ok": 0, "func_ok": 0, "exec_ok": 0, "errors": []}

    for level in config.levels:
        questions = load_questions(config.csv_path(level))
        print(f"\n  Level {level} ({len(questions)} questions)")

        for idx in range(1, len(questions) + 1):
            sol_path = config.solution_path(level, idx)
            stats["total"] += 1

            if not os.path.exists(sol_path):
                print(f"    [{idx:3d}] SKIP  MISSING")
                stats["errors"].append(f"L{level}/temp{idx}: file missing")
                continue

            # Syntax check
            syn_ok, syn_msg = check_syntax(sol_path)
            if not syn_ok:
                print(f"    [{idx:3d}] FAIL  SYNTAX: {syn_msg}")
                stats["errors"].append(f"L{level}/temp{idx}: {syn_msg}")
                continue
            stats["syntax_ok"] += 1

            # Function check
            func_ok, func_msg = check_function_exists(sol_path)
            if not func_ok:
                print(f"    [{idx:3d}] FAIL  FUNC: {func_msg}")
                stats["errors"].append(f"L{level}/temp{idx}: {func_msg}")
                continue
            stats["func_ok"] += 1

            # Quick execution test with first test input
            test_inputs = infer_test_inputs(sol_path, n_single=1)
            if test_inputs:
                args_list, kwargs = test_inputs[0]

                # Detect function name
                _, func_detail = check_function_exists(sol_path)
                func_name = "level_function"
                if "alternative:" in func_detail.lower():
                    import re
                    match = re.search(r"'(\w+)'", func_detail)
                    if match:
                        func_name = match.group(1)

                r = execute_function_direct(sol_path, func_name, args_list, kwargs, timeout=30)
                if r["success"]:
                    stats["exec_ok"] += 1
                    out_type = type(r["output"]).__name__
                    print(f"    [{idx:3d}] OK    {out_type} ({r['time_s']:.2f}s)")
                else:
                    print(f"    [{idx:3d}] FAIL  EXEC: {r['error'][:80]}")
                    stats["errors"].append(f"L{level}/temp{idx}: {r['error'][:100]}")
            else:
                print(f"    [{idx:3d}] SKIP  No test inputs inferred")

    print(f"\n{'='*60}")
    print(f"  Validation Summary")
    print(f"{'='*60}")
    print(f"  Total:     {stats['total']}")
    print(f"  Syntax OK: {stats['syntax_ok']}")
    print(f"  Func OK:   {stats['func_ok']}")
    print(f"  Exec OK:   {stats['exec_ok']}")
    print(f"  Error rate: {len(stats['errors'])}/{stats['total']}")

    if stats["errors"]:
        print(f"\n  Errors ({len(stats['errors'])}):")
        for e in stats["errors"][:20]:
            print(f"    - {e}")
        if len(stats["errors"]) > 20:
            print(f"    ... and {len(stats['errors']) - 20} more")


def cmd_generate_ground_truth(args):
    """Generate ground truth by running all solutions with standard test inputs."""
    config = EvalConfig(PROJECT_ROOT)
    config.lang = args.lang

    print("=" * 60)
    print("  Generating Ground Truth")
    print("=" * 60)

    ground_truth = {}

    for level in config.levels:
        questions = load_questions(config.csv_path(level))
        gt_level = {}

        for idx in range(1, len(questions) + 1):
            sol_path = config.solution_path(level, idx)
            if not os.path.exists(sol_path):
                continue

            test_inputs = infer_test_inputs(sol_path, config.n_test_molecules)

            # Detect function name
            _, func_detail = check_function_exists(sol_path)
            func_name = "level_function"
            if "alternative:" in func_detail.lower():
                import re
                match = re.search(r"'(\w+)'", func_detail)
                if match:
                    func_name = match.group(1)

            cases = []
            for i, (a, kw) in enumerate(test_inputs):
                r = execute_function_direct(sol_path, func_name, a, kw, timeout=30)
                cases.append({
                    "input_args": a,
                    "input_kwargs": kw,
                    "output": r["output"],
                    "success": r["success"],
                    "time_s": r["time_s"],
                })

            gt_level[idx] = {
                "question": questions[idx - 1] if idx <= len(questions) else "",
                "func_name": func_name,
                "test_cases": cases,
            }

            n_ok = sum(1 for c in cases if c["success"])
            print(f"  L{level}/temp{idx}: {n_ok}/{len(cases)} cases OK")

        ground_truth[f"level{level}"] = gt_level

    # Save
    gt_path = os.path.join(PROJECT_ROOT, "results", "ground_truth.json")
    os.makedirs(os.path.dirname(gt_path), exist_ok=True)

    def clean(obj):
        if isinstance(obj, dict):
            return {k: clean(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [clean(v) for v in obj]
        elif isinstance(obj, (int, float, str, bool, type(None))):
            return obj
        else:
            return str(obj)

    with open(gt_path, 'w', encoding='utf-8') as f:
        json.dump(clean(ground_truth), f, indent=2, ensure_ascii=False)
    print(f"\n  Ground truth saved to: {gt_path}")


def main():
    parser = argparse.ArgumentParser(
        description="RDKitBench — Molecular Vibe Coding Benchmark Evaluation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python test.py --validate-solutions         # Check all solutions
  python test.py --generate-ground-truth      # Cache ground truth
  python test.py --all                        # Evaluate all predictions
  python test.py --level 1                    # Evaluate Level 1 only
  python test.py --level 1 --question 5       # Evaluate single question
        """,
    )

    # Action
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--all", action="store_true", help="Evaluate all levels")
    group.add_argument("--level", type=int, choices=[1, 2, 3, 4, 5], help="Evaluate specific level")
    group.add_argument("--validate-solutions", action="store_true", help="Validate all solution files")
    group.add_argument("--generate-ground-truth", action="store_true", help="Generate ground truth cache")

    # Options
    parser.add_argument("--question", type=int, help="Evaluate specific question (with --level)")
    parser.add_argument("--lang", default="cn", choices=["cn", "en"], help="Question language (default: cn)")
    parser.add_argument("--timeout", type=int, default=60, help="Execution timeout in seconds")
    parser.add_argument("--n-test", type=int, default=10, help="Number of test molecules per question")
    parser.add_argument("--no-safe", action="store_true", help="Use direct execution (faster, less safe)")

    args = parser.parse_args()

    if args.all:
        cmd_evaluate_all(args)
    elif args.level:
        cmd_evaluate_level(args)
    elif args.validate_solutions:
        cmd_validate_solutions(args)
    elif args.generate_ground_truth:
        cmd_generate_ground_truth(args)


if __name__ == "__main__":
    main()
