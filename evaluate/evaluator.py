"""
RDKitBench Evaluation Framework
================================

Core evaluator: orchestrates the full evaluation pipeline.

Pipeline:
    1. Execute reference and candidate programs on shared test inputs
    2. Apply type-aware output comparison
    3. For stochastic references, try structural validation after value mismatch
    4. For executable, intrinsically non-comparable outputs, apply the strict
       AST/API-semantic fallback used by Main Pass@1
    5. Report broader API coverage separately as a diagnostic
"""

import os
import sys
import csv
import json
import time
import inspect
import re
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
from .comparators import compare_value, compare_outputs, validate_structure
from .api_matcher import api_fallback_check, extract_key_categories

# Logic test runner
try:
    from tests.runner import run_logic_test
except ImportError:
    run_logic_test = None
from .molecules import (
    SINGLE_MOLECULES,
    MOLECULE_PAIRS,
    MOLECULE_LIBRARY,
    MOLECULE_ACTIVITIES,
    MOLECULE_LABELS,
    SALT_SMILES,
    stratified_sample,
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
        self.tests_dir = os.path.join(base_dir, "tests")

        self.levels = [1, 2, 3, 4, 5]
        self.lang = "cn"  # "cn" or "en"
        self.timeout = 60
        self.float_atol = 0.01

        # Main Pass@1 uses the strict threshold only for executable outputs
        # that remain intrinsically non-comparable after Stage 2.  The broader
        # threshold is reported as a diagnostic and never changes Main Pass@1.
        self.strict_api_threshold = 0.5
        self.diagnostic_api_threshold = 0.7

        # Number of test molecules to use per question
        self.n_test_molecules = 5

        # Seed for stratified molecule sampling (different seeds → different rounds)
        self.seed = 42

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
# Evaluation-branch detection
# ============================================================

NONDETERMINISTIC_PATTERN = re.compile(
    "|".join((
        r"\brandom\.",
        r"\bnp\.random\.",
        r"\bRandomForest",
        r"\bGradientBoosting",
        r"\bKMeans\b",
        r"\bMiniBatchKMeans\b",
        r"\btrain_test_split\b",
        r"\bshuffle\b",
        r"\bsample\b",
    ))
)

NON_COMPARABLE_CATEGORIES = frozenset({
    "mol_draw_image",
    "mol_draw_svg",
    "draw_general",
    "sdf_io",
    "conformer_gen",
    "pandastools",
})

NON_COMPARABLE_DETAIL_TAGS = (
    "PIL.Image",
    "__type__",
    "fallback:",
    "ndarray",
    "rdkit.Mol",
    "SVG",
    "svg",
)


def is_nondeterministic_solution(solution_path: str) -> bool:
    """Return whether the expert solution contains a stochastic construct."""
    try:
        with open(solution_path, "r", encoding="utf-8") as handle:
            return NONDETERMINISTIC_PATTERN.search(handle.read()) is not None
    except OSError:
        return False


def is_intrinsically_non_comparable(details: list, solution_path: str) -> bool:
    """Gate the strict AST/API fallback to representation-level failures."""
    for detail in details:
        if (
            detail.get("executable")
            and not detail.get("match")
            and not detail.get("skipped")
        ):
            message = str(detail.get("detail", ""))
            if any(tag in message for tag in NON_COMPARABLE_DETAIL_TAGS):
                return True
            if "ref=" in message and any(tag in message for tag in (
                "<PIL", "<rdkit", "<class", "Image", "Mol object",
            )):
                return True

    try:
        return bool(
            extract_key_categories(solution_path) & NON_COMPARABLE_CATEGORIES
        )
    except Exception:
        return False


# ============================================================
# Test input inference
# ============================================================

def infer_test_inputs(filepath: str, n_single: int = 5, seed: int = 42) -> list:
    """
    Infer appropriate test inputs for a solution file by inspecting
    its function signature.

    For single-molecule inputs (CASE 1, the most common), uses
    stratified_sample() to draw molecules across complexity tiers.
    The seed parameter enables multiple evaluation rounds with
    different molecule samples for statistical robustness.

    Returns a list of (args, kwargs) tuples, one per test case.
    """
    # Get stratified molecules for this round
    sampled = stratified_sample(n_single, seed=seed)

    func = load_function_from_file(filepath)
    if func is None:
        # Fallback: assume single SMILES input
        return [([smi], {}) for smi in sampled]

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
        return [([smi], {}) for smi in sampled]

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
            for i in range(min(n_single, len(sampled))):
                smarts = SUBSTRUCTURE_SMARTS[i % len(SUBSTRUCTURE_SMARTS)]
                inputs.append(([sampled[i], smarts], {}))
            return inputs
        # Second molecule (e.g., scaffold hopping, MMP)
        if any(kw in second_name for kw in
               ("mol2", "smiles2", "smi2", "target", "ref")):
            return [([a, b], {}) for a, b in MOLECULE_PAIRS[:n_single]]
        # Reaction SMARTS
        if "reaction" in second_name or "rxn" in second_name:
            rxn = list(REACTION_SMARTS.values())[0]
            return [([smi, rxn], {}) for smi in sampled]
        # List of reaction SMARTS
        if "list" in second_name and "reaction" in " ".join(all_names):
            rxn_list = list(REACTION_SMARTS.values())
            return [([sampled[0], rxn_list], {})]
        # Fallback for two SMILES
        return [([a, b], {}) for a, b in MOLECULE_PAIRS[:n_single]]

    # ── CASE 3: Mol + optional seed (iterative optimization) ──
    if n_required == 1 and first_name in single_mol_names and n_total >= 2:
        # Has optional params like seed, num_confs, etc.
        return [([smi], {}) for smi in sampled]

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
                        extra_smiles = sampled
                        return [([MOLECULE_LIBRARY, MOLECULE_ACTIVITIES,
                                  extra_smiles], {})]
                return [([MOLECULE_LIBRARY, MOLECULE_ACTIVITIES], {})]
            # List + labels (QSAR classification)
            if "label" in second_name:
                if n_total >= 3:
                    third_name = all_names[2]
                    if "new" in third_name or "test" in third_name:
                        extra_smiles = sampled
                        return [([MOLECULE_LIBRARY, MOLECULE_LABELS,
                                  extra_smiles], {})]
                return [([MOLECULE_LIBRARY, MOLECULE_LABELS], {})]
            # List + target/core SMILES (BRICS recombination, R-group)
            if any(kw in second_name for kw in
                   ("target", "query", "ref", "core")):
                return [([MOLECULE_LIBRARY, sampled[0]], {})]
            # Two lists (train + test, reactants + reactants)
            if second_name in list_names or "new" in second_name:
                return [([MOLECULE_LIBRARY[:20], MOLECULE_LIBRARY[20:]], {})]
            # List + optional numeric params
            return [([MOLECULE_LIBRARY], {})]

    # ── CASE 5: Query + library (virtual screening, search) ──
    query_names = {"query_smiles", "query", "query_mol"}
    if first_name in query_names:
        if n_total >= 2:
            return [([sampled[i], MOLECULE_LIBRARY], {})
                    for i in range(min(n_single, len(sampled)))]
        return [([smi], {}) for smi in sampled]

    # ── CASE 6: Library + query (reversed order, e.g. L5/temp59) ──
    if first_name == "library_smiles" and n_total >= 2:
        second_name = all_names[1]
        if "query" in second_name:
            return [([MOLECULE_LIBRARY, sampled[i]], {})
                    for i in range(min(n_single, len(sampled)))]

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
    return [([smi], {}) for smi in sampled]


# ============================================================
# Single question evaluation
# ============================================================

def evaluate_single(solution_path: str, prediction_path: str,
                    config: EvalConfig, level: int = 0, idx: int = 0) -> dict:
    """Evaluate one prediction with the complete MolViBench cascade."""
    nondeterministic = is_nondeterministic_solution(solution_path)
    result = {
        "syntax_ok": False,
        "function_exists": False,
        "executable": False,
        "n_test_cases": 0,
        "n_pass": 0,
        "n_exact_pass": 0,
        "pass_rate": 0.0,
        "exact_pass_rate": 0.0,
        "pass_rate_api": 0.0,
        "pass_rate_fallback": 0.0,
        "main_pass": False,
        "diagnostic_fallback_pass": False,
        "nondeterministic": nondeterministic,
        "api_fallback": {"triggered": False},
        "api_fallback_broad": {"triggered": False},
        "logic_n_pass": 0,
        "logic_pass_rate": 0.0,
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
    test_inputs = infer_test_inputs(
        solution_path,
        config.n_test_molecules,
        seed=getattr(config, "seed", 42),
    )

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
    exact_pass_count = 0
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

        # Stage 2a: ordinary type-aware comparison.
        match, detail = compare_value(pred_output, ref_output, float_atol=config.float_atol)
        method = "exact"
        if match:
            pass_count += 1
            exact_pass_count += 1
        elif nondeterministic:
            # Stage 2b: stochastic references may produce different values
            # while satisfying the same output contract.  This branch does
            # not consult source-code or API overlap.
            structural_match, structural_detail = validate_structure(
                pred_output,
                ref_output,
                float_atol=config.float_atol,
            )
            if structural_match:
                match = True
                method = "structural"
                detail = structural_detail
                pass_count += 1
            else:
                method = "structural_fail"
                detail = (
                    f"exact: {detail[:100]} | "
                    f"struct: {structural_detail[:100]}"
                )

        # Run logic test
        logic_passed = None
        logic_detail = ""
        if run_logic_test is not None and level > 0 and idx > 0:
            logic_passed, logic_detail = run_logic_test(
                config.tests_dir, level, idx,
                pred_output, args, kwargs,
            )

        details.append({
            "test_case": i,
            "executable": True,
            "match": match,
            "method": method,
            "detail": detail[:200],
            "logic_passed": logic_passed,
            "logic_detail": logic_detail[:200] if logic_detail else "",
            "pred_time_s": pred_result["time_s"],
        })

    n_test = len(details)
    result["executable"] = any_executable
    result["n_test_cases"] = n_test
    result["n_pass"] = pass_count
    result["n_exact_pass"] = exact_pass_count
    result["pass_rate"] = round(pass_count / n_test, 4) if n_test > 0 else 0.0
    result["exact_pass_rate"] = (
        round(exact_pass_count / n_test, 4) if n_test > 0 else 0.0
    )

    # Logic test aggregation
    logic_tested = [d for d in details if d.get("logic_passed") is not None]
    logic_pass_count = sum(1 for d in logic_tested if d["logic_passed"])
    n_logic = len(logic_tested)
    result["logic_n_pass"] = logic_pass_count
    result["logic_n_tested"] = n_logic
    result["logic_pass_rate"] = round(logic_pass_count / n_logic, 4) if n_logic > 0 else 0.0

    result["details"] = details

    # Stage 3a (Main Pass@1): only executable, representation-level failures
    # are eligible for the strict AST/API-semantic fallback.
    non_comparable = is_intrinsically_non_comparable(details, solution_path)
    if result["pass_rate"] < 1.0 and any_executable and non_comparable:
        api_ok, api_detail = api_fallback_check(
            solution_path,
            prediction_path,
            pred_executed_ok=True,
            min_overlap_ratio=getattr(config, "strict_api_threshold", 0.5),
        )
        result["api_fallback"] = {
            "triggered": True,
            "pass": api_ok,
            "detail": api_detail,
        }
        if api_ok:
            executable_failures = sum(
                1 for item in details
                if item.get("executable")
                and not item.get("match")
                and not item.get("skipped")
            )
            result["n_pass_api"] = pass_count + executable_failures
            result["pass_rate_api"] = round(
                (pass_count + executable_failures) / n_test,
                4,
            ) if n_test else 0.0
            for item in details:
                if (
                    item.get("executable")
                    and not item.get("match")
                    and not item.get("skipped")
                ):
                    item["api_fallback_pass"] = True

    # Stage 3b: broad API coverage is a diagnostic only.  It is intentionally
    # kept separate from Main Pass@1 even when it passes.
    if result["pass_rate"] < 1.0 and any_executable:
        broad_ok, broad_detail = api_fallback_check(
            solution_path,
            prediction_path,
            pred_executed_ok=True,
            min_overlap_ratio=getattr(
                config,
                "diagnostic_api_threshold",
                0.7,
            ),
        )
        result["api_fallback_broad"] = {
            "triggered": True,
            "pass": broad_ok,
            "detail": broad_detail,
        }
        if broad_ok:
            executable_failures = sum(
                1 for item in details
                if item.get("executable")
                and not item.get("match")
                and not item.get("skipped")
            )
            result["n_pass_fallback"] = pass_count + executable_failures
            result["pass_rate_fallback"] = round(
                (pass_count + executable_failures) / n_test,
                4,
            ) if n_test else 0.0

    result["main_pass"] = (
        result["pass_rate"] >= 1.0
        or result["pass_rate_api"] >= 1.0
    )
    result["diagnostic_fallback_pass"] = (
        result["pass_rate"] >= 1.0
        or result["pass_rate_fallback"] >= 1.0
    )

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

        r = evaluate_single(sol_path, pred_path, config, level=level, idx=idx)
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
            if r.get("nondeterministic") and r.get("exact_pass_rate", 0) < 1.0:
                sym = "PASS STRUCT"
            else:
                sym = "PASS"
        elif r.get("pass_rate_api", 0) >= 1.0:
            sym = "PASS API"
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
    n_exact = sum(1 for r in evaluated if r.get("exact_pass_rate", 0) >= 1.0)
    n_main = sum(1 for r in evaluated if r.get("main_pass"))
    n_fallback = sum(
        1 for r in evaluated if r.get("diagnostic_fallback_pass")
    )
    n_partial = sum(1 for r in evaluated if 0 < r.get("pass_rate", 0) < 1.0)
    n_logic = sum(1 for r in evaluated if r.get("logic_pass_rate", 0) >= 1.0)
    n_nondeterministic = sum(
        1 for r in evaluated if r.get("nondeterministic")
    )
    n_api_triggered = sum(
        1 for r in evaluated
        if r.get("api_fallback", {}).get("triggered")
    )
    n_api_pass = sum(
        1 for r in evaluated if r.get("api_fallback", {}).get("pass")
    )
    n_broad_triggered = sum(
        1 for r in evaluated
        if r.get("api_fallback_broad", {}).get("triggered")
    )
    n_broad_pass = sum(
        1 for r in evaluated
        if r.get("api_fallback_broad", {}).get("pass")
    )

    summary = {
        "level": level,
        "n_questions": n_questions,
        "n_evaluated": n_eval,
        "n_no_prediction": n_questions - n_eval,
        "syntax_pass": n_syntax,
        "function_pass": n_func,
        "executable": n_exec,
        "exact_correct": n_exact,
        "fully_correct": n_correct,
        "fully_correct_combined": n_main,
        "main_correct": n_main,
        "fully_correct_fallback": n_fallback,
        "partially_correct": n_partial,
        "logic_correct": n_logic,
        "nondeterministic_total": n_nondeterministic,
        "api_fallback_triggered": n_api_triggered,
        "api_fallback_pass": n_api_pass,
        "broad_fallback_triggered": n_broad_triggered,
        "broad_fallback_pass": n_broad_pass,
        "rates": {
            "syntax_rate": round(n_syntax / n_eval, 4) if n_eval else 0,
            "function_rate": round(n_func / n_eval, 4) if n_eval else 0,
            "executable_rate": round(n_exec / n_eval, 4) if n_eval else 0,
            "exact_pass@1": round(n_exact / n_eval, 4) if n_eval else 0,
            "pass@1": round(n_correct / n_eval, 4) if n_eval else 0,
            "main_pass@1": round(n_main / n_eval, 4) if n_eval else 0,
            "pass@1_combined": round(n_main / n_eval, 4) if n_eval else 0,
            "fallback_pass@1": round(n_fallback / n_eval, 4) if n_eval else 0,
            "pass@1_fallback": round(n_fallback / n_eval, 4) if n_eval else 0,
            "partial_rate": round(n_partial / n_eval, 4) if n_eval else 0,
            "logic_pass@1": round(n_logic / n_eval, 4) if n_eval else 0,
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
    print("  RDKitBench Evaluation")
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
    total_exact = sum(r["exact_correct"] for r in level_results.values())
    total_correct = sum(r["fully_correct"] for r in level_results.values())
    total_main = sum(r["main_correct"] for r in level_results.values())
    total_fallback = sum(
        r["fully_correct_fallback"] for r in level_results.values()
    )
    total_logic = sum(r.get("logic_correct", 0) for r in level_results.values())
    total_api_triggered = sum(
        r.get("api_fallback_triggered", 0) for r in level_results.values()
    )
    total_api_pass = sum(
        r.get("api_fallback_pass", 0) for r in level_results.values()
    )
    total_broad_triggered = sum(
        r.get("broad_fallback_triggered", 0) for r in level_results.values()
    )
    total_broad_pass = sum(
        r.get("broad_fallback_pass", 0) for r in level_results.values()
    )

    report = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "config": {
            "lang": config.lang,
            "safe_mode": config.safe_mode,
            "n_test_molecules": config.n_test_molecules,
            "timeout": config.timeout,
            "seed": config.seed,
            "strict_api_threshold": config.strict_api_threshold,
            "diagnostic_api_threshold": config.diagnostic_api_threshold,
        },
        "grand_summary": {
            "total_questions": total_q,
            "total_evaluated": total_eval,
            "total_executable": total_exec,
            "total_exact": total_exact,
            "total_correct": total_correct,
            "total_main_correct": total_main,
            "total_correct_combined": total_main,
            "total_fallback_correct": total_fallback,
            "total_logic_correct": total_logic,
            "api_fallback_triggered": total_api_triggered,
            "api_fallback_pass": total_api_pass,
            "broad_fallback_triggered": total_broad_triggered,
            "broad_fallback_pass": total_broad_pass,
            "overall_executable_rate": round(total_exec / total_eval, 4) if total_eval else 0,
            "overall_exact_pass@1": round(total_exact / total_eval, 4) if total_eval else 0,
            "overall_pass@1": round(total_correct / total_eval, 4) if total_eval else 0,
            "overall_main_pass@1": round(total_main / total_eval, 4) if total_eval else 0,
            "overall_pass@1_combined": round(total_main / total_eval, 4) if total_eval else 0,
            "overall_fallback_pass@1": round(total_fallback / total_eval, 4) if total_eval else 0,
            "overall_pass@1_fallback": round(total_fallback / total_eval, 4) if total_eval else 0,
            "overall_logic_pass@1": round(total_logic / total_eval, 4) if total_eval else 0,
            "elapsed_s": round(elapsed, 1),
        },
        "per_level": {
            f"level{level}": {
                "n_questions": r["n_questions"],
                "executable_rate": r["rates"]["executable_rate"],
                "exact_pass@1": r["rates"]["exact_pass@1"],
                "pass@1": r["rates"]["pass@1"],
                "main_pass@1": r["rates"]["main_pass@1"],
                "fallback_pass@1": r["rates"]["fallback_pass@1"],
                "logic_pass@1": r["rates"].get("logic_pass@1", 0),
            }
            for level, r in level_results.items()
        },
        "level_details": level_results,
    }

    # Print summary table
    print(f"\n{'='*70}")
    print("  GRAND SUMMARY")
    print(f"{'='*70}")
    print(
        f"  {'Level':<8} {'Questions':<10} {'Exec':<8} {'Exact':<8} "
        f"{'Stage2':<8} {'Main':<8} {'Diag':<8}"
    )
    print(f"  {'-'*66}")
    for level, r in level_results.items():
        rates = r["rates"]
        print(f"  L{level:<7} {r['n_questions']:<10} "
              f"{rates['executable_rate']:<8.1%} "
              f"{rates['exact_pass@1']:<8.1%} "
              f"{rates['pass@1']:<8.1%} "
              f"{rates['main_pass@1']:<8.1%} "
              f"{rates['fallback_pass@1']:<8.1%}")
    print(f"  {'-'*66}")
    print(f"  {'TOTAL':<8} {total_q:<10} "
          f"{report['grand_summary']['overall_executable_rate']:<8.1%} "
          f"{report['grand_summary']['overall_exact_pass@1']:<8.1%} "
          f"{report['grand_summary']['overall_pass@1']:<8.1%} "
          f"{report['grand_summary']['overall_main_pass@1']:<8.1%} "
          f"{report['grand_summary']['overall_fallback_pass@1']:<8.1%}")
    print(
        f"\n  Strict API fallback: {total_api_pass}/{total_api_triggered} "
        "triggered -> passed"
    )
    print(
        f"  Broad diagnostic:    {total_broad_pass}/{total_broad_triggered} "
        "triggered -> passed"
    )
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
