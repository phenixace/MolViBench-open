"""
MolViBench Evaluation Framework
================================

Core evaluator: orchestrates the full evaluation pipeline.

Pipeline:
    1. Load questions from CSV
    2. Load solution code (ground truth generator)
    3. Load prediction code (LLM-generated)
    4. Run both with standard test inputs
    5. Compare outputs using type-aware comparators
    6. Aggregate metrics and generate report
"""

import os
import sys
import csv
import json
import time
import inspect
from typing import Optional
from collections import defaultdict

# Fix Windows console encoding
if sys.platform == 'win32':
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
        sys.stderr.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass

from .executor import (
    execute_function_direct,
    execute_function_subprocess,
    check_syntax,
    check_function_exists,
    load_function_from_file,
)
from .comparators import compare_value, compare_outputs
from .molecules import (
    SINGLE_MOLECULES,
    MOLECULE_PAIRS,
    MOLECULE_LIBRARY,
    MOLECULE_ACTIVITIES,
    MOLECULE_LABELS,
    SALT_SMILES,
    GENERATED_SMILES,
    SUBSTRUCTURE_SMARTS,
    RGROUP_LISTS,
    REACTION_SMARTS,
    PEPTIDE_SEQUENCES,
)


# ============================================================
# Configuration
# ============================================================

class EvalConfig:
    """Evaluation configuration."""

    def __init__(self, base_dir: str):
        self.base_dir = base_dir
        self.data_dir = os.path.join(base_dir, "data")
        self.solutions_dir = os.path.join(base_dir, "solutions")
        self.predictions_dir = os.path.join(base_dir, "predictions")
        self.results_dir = os.path.join(base_dir, "results")

        self.levels = [1, 2, 3, 4, 5]
        self.lang = "cn"  # "cn" or "en"
        self.timeout = 60
        self.float_atol = 0.01

        # Number of test molecules to use per question (from ZINC250K)
        self.n_test_molecules = 5

        # Whether to use subprocess (safe) or direct (fast) execution
        self.safe_mode = True

    def csv_path(self, level: int) -> str:
        return os.path.join(self.data_dir, self.lang, f"level{level}.csv")

    def solution_path(self, level: int, idx: int) -> str:
        return os.path.join(self.solutions_dir, f"level{level}", f"temp{idx}.py")

    def prediction_path(self, level: int, idx: int) -> str:
        return os.path.join(self.predictions_dir, f"level{level}", f"temp{idx}.py")


# ============================================================
# Question loader
# ============================================================

def load_questions(csv_path: str) -> list:
    """Load questions from a CSV file. Returns list of question strings."""
    questions = []
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            questions.append(row.get("question", "").strip())
    return questions


# ============================================================
# Test input inference
# ============================================================

def infer_test_inputs(filepath: str, n_single: int = 5) -> list:
    """
    Infer appropriate test inputs for a solution file by inspecting
    its function signature.  Uses molecules from ZINC250K.

    Returns a list of (args, kwargs) tuples, one per test case.
    Each question gets up to 5 test cases (user requirement).
    """
    func = load_function_from_file(filepath)
    if func is None:
        # Fallback: assume single SMILES input
        return [([smi], {}) for smi in SINGLE_MOLECULES[:n_single]]

    sig = inspect.signature(func)
    params = list(sig.parameters.values())

    # Filter out params with defaults only (optional params)
    required_params = [p for p in params
                       if p.default is inspect.Parameter.empty]
    all_names = [p.name.lower() for p in params]
    req_names = [p.name.lower() for p in required_params]
    n_required = len(required_params)
    n_total = len(params)

    # ── CASE 0: No parameters ──
    if n_total == 0:
        return [([], {})]

    first_name = all_names[0]

    # ── Helper: build SDF content from SMILES list ──
    def _make_sdf(smiles_list):
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            writer = Chem.SDWriter("/dev/null")  # dummy
            blocks = []
            for smi in smiles_list:
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    AllChem.Compute2DCoords(mol)
                    blocks.append(Chem.MolToMolBlock(mol))
            sdf = "\n$$$$\n".join(blocks) + "\n$$$$\n"
            return sdf
        except Exception:
            return ""

    # ── CASE 1: Single molecule (most common: ~253 files) ──
    # Names: mol, smiles, mol_smiles, smi, smiles_input, inhibitor_smiles,
    #        seed_smiles, mol_smiles
    single_mol_names = {
        "mol", "smiles", "mol_smiles", "smi", "smiles_input",
        "inhibitor_smiles", "seed_smiles", "smiles_str",
    }
    if n_required <= 1 and first_name in single_mol_names:
        return [([smi], {}) for smi in SINGLE_MOLECULES[:n_single]]

    # ── CASE 2: Two SMILES (similarity, MCS, alignment, reactions) ──
    two_mol_first = {"smiles1", "mol1", "smi1", "smiles_a",
                     "aryl_halide", "azide_smiles"}
    if first_name in two_mol_first:
        return [([a, b], {}) for a, b in MOLECULE_PAIRS[:n_single]]

    # Also catch (mol, substructure) / (mol, smarts_pattern) patterns
    if n_required == 2 and first_name in single_mol_names:
        second_name = all_names[1]
        # Substructure SMILES/SMARTS
        if any(kw in second_name for kw in
               ("sub", "smarts", "pattern", "core", "scaffold",
                "replacement", "cysteine")):
            inputs = []
            for i in range(min(n_single, len(SINGLE_MOLECULES))):
                smarts = SUBSTRUCTURE_SMARTS[i % len(SUBSTRUCTURE_SMARTS)]
                inputs.append(([SINGLE_MOLECULES[i], smarts], {}))
            return inputs
        # Second molecule (e.g., scaffold hopping, MMP)
        if any(kw in second_name for kw in
               ("mol2", "smiles2", "smi2", "target", "ref")):
            return [([a, b], {}) for a, b in MOLECULE_PAIRS[:n_single]]
        # Reaction SMARTS
        if "reaction" in second_name or "rxn" in second_name:
            rxn = list(REACTION_SMARTS.values())[0]
            return [([smi, rxn], {}) for smi in SINGLE_MOLECULES[:n_single]]
        # List of reaction SMARTS
        if "list" in second_name and "reaction" in " ".join(all_names):
            rxn_list = list(REACTION_SMARTS.values())
            return [([SINGLE_MOLECULES[0], rxn_list], {})]
        # Fallback for two SMILES
        return [([a, b], {}) for a, b in MOLECULE_PAIRS[:n_single]]

    # ── CASE 3: Mol + optional seed (iterative optimization) ──
    if n_required == 1 and first_name in single_mol_names and n_total >= 2:
        # Has optional params like seed, num_confs, etc.
        return [([smi], {}) for smi in SINGLE_MOLECULES[:n_single]]

    # ── CASE 4: SMILES list (clustering, filtering, diversity) ──
    list_names = {
        "mols", "smiles_list", "molecules", "library_smiles",
        "active_smiles_list", "train_smiles", "fragment_smiles_list",
        "library_smiles_list",
    }
    if first_name in list_names:
        if n_required == 1:
            # Just a list
            return [([MOLECULE_LIBRARY], {})]
        if n_required >= 2:
            second_name = all_names[1]
            # List + activities (QSAR regression)
            if "activit" in second_name or "pic50" in second_name:
                extra_args = {}
                # Check if there's a new_smiles_list param
                if n_total >= 3:
                    third_name = all_names[2]
                    if "new" in third_name or "test" in third_name:
                        extra_smiles = SINGLE_MOLECULES[:5]
                        return [([MOLECULE_LIBRARY, MOLECULE_ACTIVITIES,
                                  extra_smiles], {})]
                return [([MOLECULE_LIBRARY, MOLECULE_ACTIVITIES], {})]
            # List + labels (QSAR classification)
            if "label" in second_name:
                if n_total >= 3:
                    third_name = all_names[2]
                    if "new" in third_name or "test" in third_name:
                        extra_smiles = SINGLE_MOLECULES[:5]
                        return [([MOLECULE_LIBRARY, MOLECULE_LABELS,
                                  extra_smiles], {})]
                return [([MOLECULE_LIBRARY, MOLECULE_LABELS], {})]
            # List + target/core SMILES (BRICS recombination, R-group)
            if any(kw in second_name for kw in
                   ("target", "query", "ref", "core")):
                return [([MOLECULE_LIBRARY, SINGLE_MOLECULES[0]], {})]
            # Two lists (train + test, reactants + reactants)
            if second_name in list_names or "new" in second_name:
                return [([MOLECULE_LIBRARY[:20], MOLECULE_LIBRARY[20:]], {})]
            # List + optional numeric params
            return [([MOLECULE_LIBRARY], {})]

    # ── CASE 5: Query + library (virtual screening, search) ──
    query_names = {"query_smiles", "query", "query_mol"}
    if first_name in query_names:
        if n_total >= 2:
            return [([SINGLE_MOLECULES[i], MOLECULE_LIBRARY], {})
                    for i in range(min(n_single, len(SINGLE_MOLECULES)))]
        return [([smi], {}) for smi in SINGLE_MOLECULES[:n_single]]

    # ── CASE 6: Library + query (reversed order, e.g. L5/temp59) ──
    if first_name == "library_smiles" and n_total >= 2:
        second_name = all_names[1]
        if "query" in second_name:
            return [([MOLECULE_LIBRARY, SINGLE_MOLECULES[i]], {})
                    for i in range(min(n_single, len(SINGLE_MOLECULES)))]

    # ── CASE 7: Scaffold + R-groups (combinatorial) ──
    if "scaffold" in first_name:
        scaffold = "c1ccc([*:1])c([*:2])c1"
        if n_total >= 2:
            second_name = all_names[1]
            if "rgroup" in second_name or "r_group" in second_name:
                rgroup_dict = {1: RGROUP_LISTS[0], 2: RGROUP_LISTS[1]}
                return [([scaffold, rgroup_dict], {})]
            if "list" in second_name:
                return [([scaffold, RGROUP_LISTS], {})]
        return [([scaffold], {})]

    # ── CASE 8: Core SMILES + R-group dict ──
    if first_name == "core_smiles":
        rgroup_dict = {1: RGROUP_LISTS[0], 2: RGROUP_LISTS[1]}
        return [([" c1ccc([*:1])c([*:2])c1", rgroup_dict], {})]

    # ── CASE 9: Generated SMILES evaluation ──
    if "generated" in first_name:
        if n_total >= 2:
            return [([GENERATED_SMILES, MOLECULE_LIBRARY[:10]], {})]
        return [([GENERATED_SMILES], {})]

    # ── CASE 10: SMARTS pattern + library (substructure search) ──
    if first_name in ("smarts", "smarts_pattern", "pattern"):
        if n_total >= 2:
            return [([SUBSTRUCTURE_SMARTS[i], MOLECULE_LIBRARY], {})
                    for i in range(min(n_single, len(SUBSTRUCTURE_SMARTS)))]
        return [([s], {}) for s in SUBSTRUCTURE_SMARTS[:n_single]]

    # ── CASE 11: Reactants lists (reaction tasks) ──
    if first_name == "reactants" or first_name == "reactant_smiles":
        if n_total >= 2:
            second_name = all_names[1]
            if "reaction" in second_name or "rxn" in second_name:
                rxn = list(REACTION_SMARTS.values())[0]
                return [([SINGLE_MOLECULES[:3], rxn], {})]
            if "product" in second_name:
                return [([SINGLE_MOLECULES[:3], SINGLE_MOLECULES[3:6]], {})]
        return [([SINGLE_MOLECULES[:3]], {})]

    # ── CASE 12: Reactant lists for enumerate_products (L2/temp69) ──
    if first_name in ("reactants_list1", "reactant_list_1"):
        return [([SINGLE_MOLECULES[:3], SINGLE_MOLECULES[3:6],
                  list(REACTION_SMARTS.values())[0]], {})]

    # ── CASE 13: Reaction SMILES string (e.g. "A.B>>C") ──
    if first_name in ("reaction_smiles", "rxn_smiles", "reaction"):
        rxn_smi = f"{SINGLE_MOLECULES[0]}.{SINGLE_MOLECULES[1]}>>{SINGLE_MOLECULES[2]}"
        return [([rxn_smi], {})]

    # ── CASE 14: Peptide / amino acid sequence ──
    if first_name in ("sequence", "aa_sequence", "peptide", "peptide_sequence"):
        return [([seq], {}) for seq in PEPTIDE_SEQUENCES[:n_single]]

    # ── CASE 15: SDF content string ──
    if first_name in ("sdf_content", "sdf_string", "sdf_data", "sdf_text"):
        sdf = _make_sdf(MOLECULE_LIBRARY[:10])
        return [([sdf], {})]

    # ── CASE 16: SMILES list (as mols_smiles, etc.) ──
    if "list" in first_name or "mols" in first_name:
        return [([MOLECULE_LIBRARY], {})]

    # ── CASE 17: Target MW range (L4/temp57 — 0 required params) ──
    if n_required == 0:
        return [([], {})]

    # ── DEFAULT: assume single SMILES input ──
    return [([smi], {}) for smi in SINGLE_MOLECULES[:n_single]]


# ============================================================
# Single question evaluation
# ============================================================

def evaluate_single(solution_path: str, prediction_path: str,
                    config: EvalConfig) -> dict:
    """
    Evaluate a single prediction against its solution.

    Returns:
        {
            "question_idx": int,
            "syntax_ok": bool,
            "function_exists": bool,
            "executable": bool,
            "n_test_cases": int,
            "n_pass": int,
            "pass_rate": float,
            "errors": list,
            "details": list,
        }
    """
    result = {
        "syntax_ok": False,
        "function_exists": False,
        "executable": False,
        "n_test_cases": 0,
        "n_pass": 0,
        "pass_rate": 0.0,
        "errors": [],
        "details": [],
    }

    # Check if prediction file exists
    if not os.path.exists(prediction_path):
        result["errors"].append("Prediction file not found")
        return result

    # Step 1: Syntax check
    syn_ok, syn_msg = check_syntax(prediction_path)
    result["syntax_ok"] = syn_ok
    if not syn_ok:
        result["errors"].append(f"Syntax error: {syn_msg}")
        return result

    # Step 2: Function existence check
    func_ok, func_msg = check_function_exists(prediction_path)
    result["function_exists"] = func_ok
    if not func_ok:
        result["errors"].append(f"Function missing: {func_msg}")
        return result

    # Step 3: Determine test inputs from solution
    test_inputs = infer_test_inputs(solution_path, config.n_test_molecules)

    # Step 4: Run solution to get ground truth
    exec_fn = execute_function_subprocess if config.safe_mode else execute_function_direct

    # Detect the function name in prediction
    _, func_detail = check_function_exists(prediction_path)
    pred_func_name = "level_function"  # default
    if "alternative:" in func_detail.lower():
        # Extract the actual function name found
        import re
        match = re.search(r"'(\w+)'", func_detail)
        if match:
            pred_func_name = match.group(1)

    # Also detect solution function name
    _, sol_func_detail = check_function_exists(solution_path)
    sol_func_name = "level_function"
    if "alternative:" in sol_func_detail.lower():
        import re
        match = re.search(r"'(\w+)'", sol_func_detail)
        if match:
            sol_func_name = match.group(1)

    any_executable = False
    pass_count = 0
    details = []

    for i, (args, kwargs) in enumerate(test_inputs):
        # Run solution (trusted — direct execution is OK)
        sol_result = execute_function_direct(solution_path, sol_func_name, args, kwargs, config.timeout)

        if not sol_result["success"]:
            details.append({
                "test_case": i,
                "skipped": True,
                "reason": f"Solution failed: {sol_result['error']}",
            })
            continue

        ref_output = sol_result["output"]

        # Run prediction (untrusted — subprocess)
        pred_result = exec_fn(prediction_path, pred_func_name, args, kwargs, config.timeout)

        if not pred_result["success"]:
            details.append({
                "test_case": i,
                "executable": False,
                "match": False,
                "error": pred_result.get("error", "Unknown error")[:200],
            })
            continue

        any_executable = True
        pred_output = pred_result["output"]

        # Compare outputs
        match, detail = compare_value(pred_output, ref_output, float_atol=config.float_atol)
        if match:
            pass_count += 1

        details.append({
            "test_case": i,
            "executable": True,
            "match": match,
            "detail": detail[:200],
            "pred_time_s": pred_result["time_s"],
        })

    n_test = len(details)
    result["executable"] = any_executable
    result["n_test_cases"] = n_test
    result["n_pass"] = pass_count
    result["pass_rate"] = round(pass_count / n_test, 4) if n_test > 0 else 0.0
    result["details"] = details

    return result


# ============================================================
# Level evaluation
# ============================================================

def evaluate_level(level: int, config: EvalConfig) -> dict:
    """Evaluate all questions in a given level."""
    questions = load_questions(config.csv_path(level))
    n_questions = len(questions)

    print(f"\n{'='*60}")
    print(f"  Level {level}: {n_questions} questions")
    print(f"{'='*60}")

    results = []
    for idx in range(1, n_questions + 1):
        sol_path = config.solution_path(level, idx)
        pred_path = config.prediction_path(level, idx)

        print(f"  [{idx:3d}/{n_questions}] ", end="", flush=True)

        if not os.path.exists(sol_path):
            print(f"SKIP (no solution)")
            results.append({"question_idx": idx, "skipped": True})
            continue

        r = evaluate_single(sol_path, pred_path, config)
        r["question_idx"] = idx
        r["question"] = questions[idx - 1] if idx <= len(questions) else ""
        results.append(r)

        # Status symbol
        if r.get("skipped"):
            sym = "SKIP"
        elif not r["syntax_ok"]:
            sym = "FAIL SYN"
        elif not r["function_exists"]:
            sym = "FAIL FUN"
        elif not r["executable"]:
            sym = "FAIL EXE"
        elif r["pass_rate"] >= 1.0:
            sym = "PASS"
        elif r["pass_rate"] > 0:
            sym = f"PART {r['pass_rate']:.0%}"
        else:
            sym = "FAIL OUT"

        print(f"{sym}  {r.get('question', '')[:50]}")

    # Aggregate
    evaluated = [r for r in results if not r.get("skipped")]
    n_eval = len(evaluated)
    n_syntax = sum(1 for r in evaluated if r.get("syntax_ok"))
    n_func = sum(1 for r in evaluated if r.get("function_exists"))
    n_exec = sum(1 for r in evaluated if r.get("executable"))
    n_correct = sum(1 for r in evaluated if r.get("pass_rate", 0) >= 1.0)
    n_partial = sum(1 for r in evaluated if 0 < r.get("pass_rate", 0) < 1.0)

    summary = {
        "level": level,
        "n_questions": n_questions,
        "n_evaluated": n_eval,
        "n_no_prediction": n_questions - n_eval,
        "syntax_pass": n_syntax,
        "function_pass": n_func,
        "executable": n_exec,
        "fully_correct": n_correct,
        "partially_correct": n_partial,
        "rates": {
            "syntax_rate": round(n_syntax / n_eval, 4) if n_eval else 0,
            "function_rate": round(n_func / n_eval, 4) if n_eval else 0,
            "executable_rate": round(n_exec / n_eval, 4) if n_eval else 0,
            "pass@1": round(n_correct / n_eval, 4) if n_eval else 0,
            "partial_rate": round(n_partial / n_eval, 4) if n_eval else 0,
        },
        "per_question": results,
    }

    return summary


# ============================================================
# Full benchmark evaluation
# ============================================================

def evaluate_all(config: EvalConfig) -> dict:
    """Run evaluation across all levels."""
    print("=" * 60)
    print("  MolViBench Evaluation")
    print(f"  Language: {config.lang} | Safe mode: {config.safe_mode}")
    print(f"  Test molecules per question: {config.n_test_molecules}")
    print("=" * 60)

    start_time = time.time()
    level_results = {}

    for level in config.levels:
        level_results[level] = evaluate_level(level, config)

    elapsed = time.time() - start_time

    # Grand summary
    total_q = sum(r["n_questions"] for r in level_results.values())
    total_eval = sum(r["n_evaluated"] for r in level_results.values())
    total_exec = sum(r["executable"] for r in level_results.values())
    total_correct = sum(r["fully_correct"] for r in level_results.values())

    report = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "config": {
            "lang": config.lang,
            "safe_mode": config.safe_mode,
            "n_test_molecules": config.n_test_molecules,
            "timeout": config.timeout,
        },
        "grand_summary": {
            "total_questions": total_q,
            "total_evaluated": total_eval,
            "total_executable": total_exec,
            "total_correct": total_correct,
            "overall_executable_rate": round(total_exec / total_eval, 4) if total_eval else 0,
            "overall_pass@1": round(total_correct / total_eval, 4) if total_eval else 0,
            "elapsed_s": round(elapsed, 1),
        },
        "per_level": {
            f"level{level}": {
                "n_questions": r["n_questions"],
                "executable_rate": r["rates"]["executable_rate"],
                "pass@1": r["rates"]["pass@1"],
            }
            for level, r in level_results.items()
        },
        "level_details": level_results,
    }

    # Print summary table
    print(f"\n{'='*60}")
    print("  GRAND SUMMARY")
    print(f"{'='*60}")
    print(f"  {'Level':<8} {'Questions':<10} {'Syntax':<8} {'Func':<8} {'Exec':<8} {'Pass@1':<8}")
    print(f"  {'-'*52}")
    for level, r in level_results.items():
        rates = r["rates"]
        print(f"  L{level:<7} {r['n_questions']:<10} "
              f"{rates['syntax_rate']:<8.1%} "
              f"{rates['function_rate']:<8.1%} "
              f"{rates['executable_rate']:<8.1%} "
              f"{rates['pass@1']:<8.1%}")
    print(f"  {'-'*52}")
    print(f"  {'TOTAL':<8} {total_q:<10} "
          f"{'':8} {'':8} "
          f"{report['grand_summary']['overall_executable_rate']:<8.1%} "
          f"{report['grand_summary']['overall_pass@1']:<8.1%}")
    print(f"\n  Time: {elapsed:.1f}s")

    return report


# ============================================================
# Save / Load results
# ============================================================

def save_report(report: dict, filepath: str):
    """Save evaluation report as JSON."""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    # Remove non-serializable parts
    def clean(obj):
        if isinstance(obj, dict):
            return {k: clean(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [clean(v) for v in obj]
        elif isinstance(obj, (int, float, str, bool, type(None))):
            return obj
        else:
            return str(obj)

    with open(filepath, 'w', encoding='utf-8') as f:
        json.dump(clean(report), f, indent=2, ensure_ascii=False)
    print(f"\n  Report saved to: {filepath}")
