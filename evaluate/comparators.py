"""
RDKitBench Evaluation Framework
================================

Type-aware output comparators for deterministic evaluation.
Each comparator returns (is_match: bool, detail: str).
"""

from rdkit import Chem
import math
import numpy as np


# ============================================================
# Atomic comparators
# ============================================================

def compare_int(pred, ref):
    """Exact integer match."""
    try:
        p, r = int(pred), int(ref)
        match = (p == r)
        return match, f"pred={p}, ref={r}"
    except (TypeError, ValueError) as e:
        return False, f"Type error: {e}"


def compare_float(pred, ref, atol=1e-2, rtol=1e-3):
    """Numeric match within absolute + relative tolerance."""
    try:
        p, r = float(pred), float(ref)
        if math.isnan(p) and math.isnan(r):
            return True, "both NaN"
        if math.isnan(p) or math.isnan(r):
            return False, f"pred={p}, ref={r} (one is NaN)"
        abs_ok = abs(p - r) <= atol
        rel_ok = abs(p - r) <= rtol * max(abs(r), 1e-10)
        match = abs_ok or rel_ok
        return match, f"pred={p:.6f}, ref={r:.6f}, diff={abs(p-r):.2e}"
    except (TypeError, ValueError) as e:
        return False, f"Type error: {e}"


def compare_bool(pred, ref):
    """Exact boolean match."""
    try:
        p, r = bool(pred), bool(ref)
        match = (p == r)
        return match, f"pred={p}, ref={r}"
    except (TypeError, ValueError) as e:
        return False, f"Type error: {e}"


def compare_smiles(pred, ref):
    """Canonical SMILES equivalence."""
    try:
        if pred is None or ref is None:
            match = (pred is None and ref is None)
            return match, f"pred={pred}, ref={ref}"

        mol_p = Chem.MolFromSmiles(str(pred))
        mol_r = Chem.MolFromSmiles(str(ref))

        if mol_p is None and mol_r is None:
            return True, "both invalid SMILES"
        if mol_p is None or mol_r is None:
            return False, f"pred_valid={mol_p is not None}, ref_valid={mol_r is not None}"

        can_p = Chem.MolToSmiles(mol_p)
        can_r = Chem.MolToSmiles(mol_r)
        match = (can_p == can_r)
        return match, f"pred_can={can_p}, ref_can={can_r}"
    except Exception as e:
        return False, f"Error: {e}"


def compare_str(pred, ref):
    """Exact string match (after strip). Falls back to SMILES comparison if both parseable."""
    try:
        p, r = str(pred).strip(), str(ref).strip()
        if p == r:
            return True, f"exact match"

        # Try SMILES equivalence as fallback
        mol_p = Chem.MolFromSmiles(p)
        mol_r = Chem.MolFromSmiles(r)
        if mol_p is not None and mol_r is not None:
            can_p = Chem.MolToSmiles(mol_p)
            can_r = Chem.MolToSmiles(mol_r)
            if can_p == can_r:
                return True, f"SMILES equivalent: {can_p}"

        return False, f"pred='{p[:80]}', ref='{r[:80]}'"
    except Exception as e:
        return False, f"Error: {e}"


# ============================================================
# Collection comparators
# ============================================================

def compare_list(pred, ref, element_comparator=None):
    """Sorted-set comparison for lists. Handles list[str], list[int], list[float]."""
    try:
        if not isinstance(pred, (list, tuple)) or not isinstance(ref, (list, tuple)):
            return False, f"Type mismatch: pred={type(pred).__name__}, ref={type(ref).__name__}"

        if len(pred) != len(ref):
            return False, f"Length mismatch: pred={len(pred)}, ref={len(ref)}"

        if len(pred) == 0:
            return True, "both empty"

        # Determine element type from reference
        sample = ref[0]

        if isinstance(sample, str):
            # Try SMILES canonical sort
            def canonicalize(s):
                mol = Chem.MolFromSmiles(s)
                return Chem.MolToSmiles(mol) if mol else s

            try:
                pred_sorted = sorted(canonicalize(str(x)) for x in pred)
                ref_sorted = sorted(canonicalize(str(x)) for x in ref)
            except Exception:
                pred_sorted = sorted(str(x) for x in pred)
                ref_sorted = sorted(str(x) for x in ref)

            match = (pred_sorted == ref_sorted)
            n_match = sum(1 for a, b in zip(pred_sorted, ref_sorted) if a == b)
            return match, f"{n_match}/{len(ref)} elements match"

        elif isinstance(sample, (int, bool)):
            pred_sorted = sorted(int(x) for x in pred)
            ref_sorted = sorted(int(x) for x in ref)
            match = (pred_sorted == ref_sorted)
            return match, f"match={match}"

        elif isinstance(sample, float):
            pred_sorted = sorted(float(x) for x in pred)
            ref_sorted = sorted(float(x) for x in ref)
            diffs = [abs(a - b) for a, b in zip(pred_sorted, ref_sorted)]
            max_diff = max(diffs)
            match = (max_diff < 0.01)
            return match, f"max_diff={max_diff:.4e}"

        else:
            # Fallback: string comparison
            pred_sorted = sorted(str(x) for x in pred)
            ref_sorted = sorted(str(x) for x in ref)
            match = (pred_sorted == ref_sorted)
            return match, f"match={match}"

    except Exception as e:
        return False, f"Error: {e}"


def compare_dict(pred, ref, float_atol=0.01):
    """Recursive key-value comparison for dicts."""
    try:
        if not isinstance(pred, dict) or not isinstance(ref, dict):
            return False, f"Type mismatch: pred={type(pred).__name__}, ref={type(ref).__name__}"

        # Check key sets
        pred_keys = set(pred.keys())
        ref_keys = set(ref.keys())
        missing = ref_keys - pred_keys
        extra = pred_keys - ref_keys

        if missing:
            return False, f"Missing keys: {missing}"

        # Compare each value in ref
        mismatches = []
        for key in ref_keys:
            p_val, r_val = pred.get(key), ref[key]
            ok, detail = compare_value(p_val, r_val, float_atol=float_atol)
            if not ok:
                mismatches.append(f"{key}: {detail}")

        if mismatches:
            return False, f"{len(mismatches)} mismatches: {'; '.join(mismatches[:5])}"

        return True, f"all {len(ref_keys)} keys match"

    except Exception as e:
        return False, f"Error: {e}"


# ============================================================
# Universal dispatcher
# ============================================================

def compare_value(pred, ref, float_atol=0.01):
    """Auto-detect type and dispatch to appropriate comparator."""
    # Handle None
    if ref is None:
        return (pred is None), f"pred={pred}, ref=None"
    if pred is None:
        return False, f"pred=None, ref={ref}"

    # ── Handle smart_serialize'd objects (from subprocess) ──
    # If pred is a dict with __type__, it was a non-JSON-serializable object
    # serialized by the subprocess wrapper. Try type-aware comparison.
    if isinstance(pred, dict) and "__type__" in pred:
        return _compare_typed_object(pred, ref, float_atol)
    # If ref is a typed object (solution ran direct, returned special object)
    if isinstance(ref, dict) and "__type__" in ref:
        return _compare_typed_object(pred, ref, float_atol)

    # Type-based dispatch (use reference type as ground truth)
    if isinstance(ref, bool):
        return compare_bool(pred, ref)
    elif isinstance(ref, int) and not isinstance(ref, bool):
        # Allow float that equals int
        if isinstance(pred, float) and pred == int(pred):
            return compare_int(int(pred), ref)
        return compare_int(pred, ref)
    elif isinstance(ref, float):
        return compare_float(pred, ref, atol=float_atol)
    elif isinstance(ref, str):
        return compare_str(pred, ref)
    elif isinstance(ref, set):
        # Convert sets to sorted lists for comparison
        return compare_list(sorted(str(x) for x in pred) if isinstance(pred, (set, list, tuple)) else pred,
                           sorted(str(x) for x in ref))
    elif isinstance(ref, dict):
        return compare_dict(pred, ref, float_atol=float_atol)
    elif isinstance(ref, (list, tuple)):
        return compare_list(pred, ref)
    else:
        # Check if ref is a non-serializable object (PIL.Image, Mol, etc.)
        ref_type = type(ref).__name__
        ref_module = type(ref).__module__ or ""

        # PIL Image: if pred is also an image-like thing, consider it a structural match
        if "PIL" in ref_module or ref_type == "Image":
            if isinstance(pred, dict) and pred.get("__type__") == "PIL.Image":
                return True, f"both PIL.Image (pred serialized: {pred})"
            # pred might also be a PIL Image (both ran direct)
            if "PIL" in str(type(pred).__module__):
                return True, "both PIL.Image objects"
            return False, f"ref=PIL.Image, pred={type(pred).__name__}"

        # numpy array
        if ref_type == "ndarray":
            try:
                import numpy as np
                if isinstance(pred, (list, tuple)):
                    pred_arr = np.array(pred)
                elif isinstance(pred, dict) and "values" in pred:
                    pred_arr = np.array(pred["values"])
                elif hasattr(pred, 'shape'):
                    pred_arr = pred
                else:
                    return False, f"ref=ndarray, pred={type(pred).__name__}"
                return bool(np.allclose(pred_arr, ref, atol=float_atol)), \
                    f"ndarray shape ref={ref.shape}, pred={pred_arr.shape}"
            except Exception as e:
                return False, f"ndarray compare error: {e}"

        # Fallback to string comparison
        return str(pred) == str(ref), f"fallback: pred={str(pred)[:50]}, ref={str(ref)[:50]}"


def _compare_typed_object(pred, ref, float_atol=0.01):
    """Compare when one or both sides are smart_serialize'd typed objects."""
    pred_type = pred.get("__type__", "") if isinstance(pred, dict) else ""
    ref_type_name = type(ref).__name__
    ref_module = type(ref).__module__ or ""

    # Both are typed dicts (both from subprocess)
    if isinstance(pred, dict) and "__type__" in pred and isinstance(ref, dict) and "__type__" in ref:
        if pred["__type__"] == ref["__type__"]:
            # Same type — compare available fields
            if "smiles" in pred and "smiles" in ref:
                return compare_smiles(pred["smiles"], ref["smiles"])
            if "values" in pred and "values" in ref:
                return compare_list(pred["values"], ref["values"])
            return True, f"same typed object: {pred['__type__']}"
        return False, f"type mismatch: pred={pred['__type__']}, ref={ref['__type__']}"

    # pred is typed dict, ref is raw object
    if isinstance(pred, dict) and "__type__" in pred:
        # PIL Image
        if pred_type == "PIL.Image" and ("PIL" in ref_module or ref_type_name == "Image"):
            return True, f"both PIL.Image (pred serialized)"
        # numpy array
        if pred_type == "numpy.ndarray" and ref_type_name == "ndarray":
            if "values" in pred:
                try:
                    import numpy as np
                    return bool(np.allclose(np.array(pred["values"]), ref, atol=float_atol)), \
                        f"ndarray compare"
                except Exception:
                    pass
            return True, f"both ndarray (structural match)"
        # RDKit Mol
        if pred_type == "rdkit.Mol" and hasattr(ref, "GetNumAtoms"):
            if "smiles" in pred:
                try:
                    ref_smi = Chem.MolToSmiles(ref)
                    return compare_smiles(pred["smiles"], ref_smi)
                except Exception:
                    pass
            return True, f"both rdkit.Mol"
        # Generic: type name match is enough for structural equivalence
        if ref_type_name in pred_type or pred_type in str(type(ref)):
            return True, f"type match: {pred_type}"
        return False, f"typed mismatch: pred={pred_type}, ref={ref_type_name}"

    # ref is typed dict, pred is raw — symmetric case
    if isinstance(ref, dict) and "__type__" in ref:
        return _compare_typed_object(ref, pred, float_atol)

    return False, f"unexpected typed comparison"


# ============================================================
# Structural validation for stochastic reference solutions
# ============================================================

def validate_structure(pred, ref, float_atol=0.01):
    """Validate the output contract without requiring identical values.

    This comparator is used only after ordinary type-aware comparison fails
    for a task whose expert solution is flagged as stochastic.  It checks the
    returned container/schema and basic molecular validity; it is deliberately
    separate from the AST/API fallback.
    """
    if ref is None:
        return pred is None, "ref=None"
    if pred is None:
        return False, "pred=None, ref is not None"

    ref_type = type(ref).__name__
    pred_type = type(pred).__name__

    if isinstance(ref, dict):
        if not isinstance(pred, dict):
            return False, f"Type mismatch: pred={pred_type}, ref=dict"

        ref_keys = set(ref)
        pred_keys = set(pred)
        shared_keys = ref_keys & pred_keys
        if ref_keys:
            overlap = len(shared_keys) / len(ref_keys)
            if overlap < 0.5:
                return False, (
                    f"Key overlap too low: {overlap:.0%} ({shared_keys})"
                )

        mismatches = []
        for key in shared_keys:
            ok, detail = validate_structure(pred[key], ref[key], float_atol)
            if not ok:
                mismatches.append(f"{key}: {detail}")
        if len(mismatches) > len(ref_keys) * 0.5:
            return False, f"Too many structural mismatches: {mismatches[:3]}"
        return True, f"dict structure OK ({len(ref_keys)} keys)"

    if isinstance(ref, (list, tuple)):
        if not isinstance(pred, (list, tuple)):
            return False, f"Type mismatch: pred={pred_type}, ref=list"
        if not ref:
            return len(pred) == 0, "empty list check"
        if not pred:
            return False, "pred list is empty, ref is not"
        ok, detail = validate_structure(pred[0], ref[0], float_atol)
        return ok, f"list[{len(pred)}] first elem: {detail}"

    if isinstance(ref, bool):
        return isinstance(pred, (bool, int)), f"bool check: pred={pred}"

    if isinstance(ref, (int, float)):
        try:
            float(pred)
            return True, f"numeric: pred={pred}"
        except (TypeError, ValueError):
            return False, f"pred not numeric: {pred}"

    if isinstance(ref, str):
        if not isinstance(pred, str):
            return False, f"Type mismatch: pred={pred_type}, ref=str"
        if pred:
            pred_mol = Chem.MolFromSmiles(pred)
            ref_mol = Chem.MolFromSmiles(ref)
            if ref_mol is not None:
                return (
                    pred_mol is not None,
                    f"SMILES validity: pred_valid={pred_mol is not None}",
                )
        return True, "string OK"

    # For uncommon opaque values, preserve the original structural policy:
    # successful execution plus a compatible outer contract is sufficient.
    return True, f"type={ref_type} (fallback OK)"


# ============================================================
# Aggregate comparison for a test case
# ============================================================

def compare_outputs(pred_outputs, ref_outputs, float_atol=0.01):
    """
    Compare a list of (pred, ref) output pairs.
    Returns overall pass rate and per-case details.
    """
    results = []
    for i, (pred, ref) in enumerate(zip(pred_outputs, ref_outputs)):
        match, detail = compare_value(pred, ref, float_atol=float_atol)
        results.append({
            "test_case": i,
            "match": match,
            "detail": detail,
            "pred_type": type(pred).__name__,
            "ref_type": type(ref).__name__,
        })

    n_pass = sum(1 for r in results if r["match"])
    n_total = len(results)
    pass_rate = n_pass / n_total if n_total > 0 else 0.0

    return {
        "pass_rate": round(pass_rate, 4),
        "n_pass": n_pass,
        "n_total": n_total,
        "details": results,
    }
