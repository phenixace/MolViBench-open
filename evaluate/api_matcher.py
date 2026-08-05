"""
MolViBench Evaluation Framework
================================

API-based fallback evaluator.

When output comparison fails, this module checks whether the prediction
code calls the same key RDKit/cheminformatics APIs as the solution.

Logic:
    1. Extract "key API calls" from solution code via AST analysis
    2. Extract the same from prediction code
    3. If the prediction calls a sufficient subset of the solution's
       key APIs AND executes without error, consider it a pass.

This handles cases where:
    - The prediction returns a different type (Image vs SVG string)
    - The prediction wraps the result differently (dict vs raw value)
    - The prediction does the right work but returns None / prints instead
    - The return value is non-deterministic (3D coords, MolBlock)
"""

import ast
import os
from typing import Set, Tuple, Dict, List, Optional


# ============================================================
# RDKit API taxonomy — grouped by functional category
# ============================================================
# Each key is a category name; value is a set of qualified names
# that indicate the code is performing that category of work.
# We match at the attribute level (e.g. "Chem.MolFromSmiles").

RDKIT_API_CATEGORIES = {
    # ── Molecule I/O ──
    "mol_from_smiles": {
        "Chem.MolFromSmiles", "MolFromSmiles",
    },
    "mol_from_smarts": {
        "Chem.MolFromSmarts", "MolFromSmarts",
    },
    "mol_to_smiles": {
        "Chem.MolToSmiles", "MolToSmiles",
    },
    "mol_to_smarts": {
        "Chem.MolToSmarts", "MolToSmarts",
    },
    "mol_from_molblock": {
        "Chem.MolFromMolBlock", "MolFromMolBlock",
    },
    "mol_to_molblock": {
        "Chem.MolToMolBlock", "MolToMolBlock",
    },
    "mol_from_pdb": {
        "Chem.MolFromPDBBlock", "Chem.MolFromPDBFile", "MolFromPDBBlock", "MolFromPDBFile",
    },
    "mol_to_pdb": {
        "Chem.MolToPDBBlock", "Chem.MolToPDBFile", "MolToPDBBlock", "MolToPDBFile",
    },
    "mol_from_mol2": {
        "Chem.MolFromMol2Block", "Chem.MolFromMol2File",
    },
    "sdf_io": {
        "Chem.SDWriter", "SDWriter", "Chem.SDMolSupplier", "SDMolSupplier",
        "Chem.ForwardSDMolSupplier", "ForwardSDMolSupplier",
    },

    # ── 2D/3D Visualization ──
    "mol_draw_image": {
        "Draw.MolToImage", "MolToImage",
        "Draw.MolsToGridImage", "MolsToGridImage",
    },
    "mol_draw_svg": {
        "rdMolDraw2D.MolDraw2DSVG", "MolDraw2DSVG",
        "rdMolDraw2D.MolDraw2DCairo", "MolDraw2DCairo",
        "Draw.MolToFile", "MolToFile",
    },
    "draw_general": {
        "DrawMolecule", "drawer.DrawMolecule",
        "FinishDrawing", "drawer.FinishDrawing",
        "GetDrawingText", "drawer.GetDrawingText",
    },

    # ── Descriptors ──
    "descriptors_basic": {
        "Descriptors.MolWt", "Descriptors.ExactMolWt",
        "Descriptors.MolLogP", "Descriptors.TPSA",
        "Descriptors.NumHDonors", "Descriptors.NumHAcceptors",
        "Descriptors.NumRotatableBonds",
        "Descriptors.NumHeteroatoms",
        "Descriptors.RingCount",
        "Descriptors.FractionCSP3",
        "Descriptors.HeavyAtomCount",
        "Descriptors.NumAromaticRings",
        "Descriptors.NumAliphaticRings",
        "Descriptors.LabuteASA",
        "Descriptors.BertzCT",
        "rdMolDescriptors.CalcExactMolWt", "CalcExactMolWt",
        "rdMolDescriptors.CalcMolFormula", "CalcMolFormula",
        "rdMolDescriptors.CalcTPSA", "CalcTPSA",
        "rdMolDescriptors.CalcNumRotatableBonds", "CalcNumRotatableBonds",
        "rdMolDescriptors.CalcNumHBA", "CalcNumHBA",
        "rdMolDescriptors.CalcNumHBD", "CalcNumHBD",
        "rdMolDescriptors.CalcNumRings", "CalcNumRings",
        "rdMolDescriptors.CalcNumAromaticRings", "CalcNumAromaticRings",
        "rdMolDescriptors.CalcFractionCSP3", "CalcFractionCSP3",
        "rdMolDescriptors.CalcNumHeteroatoms", "CalcNumHeteroatoms",
        "rdMolDescriptors.CalcLabuteASA", "CalcLabuteASA",
        "rdMolDescriptors.CalcNumAmideBonds", "CalcNumAmideBonds",
        "rdMolDescriptors.CalcNumSpiroAtoms", "CalcNumSpiroAtoms",
        "rdMolDescriptors.CalcNumBridgeheadAtoms", "CalcNumBridgeheadAtoms",
        "rdMolDescriptors.CalcNumAtomStereoCenters", "CalcNumAtomStereoCenters",
        "rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters",
    },
    "descriptors_lipinski": {
        "Lipinski", "Descriptors.NumHDonors", "Descriptors.NumHAcceptors",
        "Descriptors.MolWt", "Descriptors.MolLogP",
        "FilterCatalog",
    },
    "descriptors_qed": {
        "QED.qed", "QED.default", "QED.properties",
    },

    # ── Fingerprints ──
    "fp_morgan": {
        "AllChem.GetMorganFingerprintAsBitVect",
        "AllChem.GetMorganFingerprint",
        "rdMolDescriptors.GetMorganFingerprintAsBitVect",
        "GetMorganFingerprintAsBitVect",
        "GetMorganFingerprint",
    },
    "fp_rdkit": {
        "Chem.RDKFingerprint", "RDKFingerprint",
        "rdmolops.RDKFingerprint",
    },
    "fp_maccs": {
        "MACCSkeys.GenMACCSKeys", "GenMACCSKeys",
    },
    "fp_topological": {
        "rdMolDescriptors.GetAtomPairFingerprint",
        "rdMolDescriptors.GetTopologicalTorsionFingerprint",
        "GetAtomPairFingerprint", "GetTopologicalTorsionFingerprint",
    },

    # ── Similarity ──
    "similarity": {
        "DataStructs.TanimotoSimilarity", "TanimotoSimilarity",
        "DataStructs.DiceSimilarity", "DiceSimilarity",
        "DataStructs.BulkTanimotoSimilarity", "BulkTanimotoSimilarity",
        "DataStructs.FingerprintSimilarity", "FingerprintSimilarity",
    },

    # ── Substructure ──
    "substructure": {
        "HasSubstructMatch", "mol.HasSubstructMatch",
        "GetSubstructMatch", "mol.GetSubstructMatch",
        "GetSubstructMatches", "mol.GetSubstructMatches",
    },

    # ── 3D Conformer ──
    "conformer_gen": {
        "AllChem.EmbedMolecule", "EmbedMolecule",
        "AllChem.EmbedMultipleConfs", "EmbedMultipleConfs",
        "AllChem.ETKDGv3", "ETKDGv3", "AllChem.ETKDG", "ETKDG",
        "rdDistGeom.EmbedMolecule", "rdDistGeom.EmbedMultipleConfs",
    },
    "conformer_optimize": {
        "AllChem.MMFFOptimizeMolecule", "MMFFOptimizeMolecule",
        "AllChem.UFFOptimizeMolecule", "UFFOptimizeMolecule",
        "AllChem.MMFFOptimizeMoleculeConfs", "MMFFOptimizeMoleculeConfs",
        "AllChem.UFFOptimizeMoleculeConfs", "UFFOptimizeMoleculeConfs",
    },
    "conformer_energy": {
        "AllChem.MMFFGetMoleculeForceField", "MMFFGetMoleculeForceField",
        "AllChem.UFFGetMoleculeForceField", "UFFGetMoleculeForceField",
        "ForceField.CalcEnergy", "CalcEnergy",
        "AllChem.MMFFGetMoleculeProperties", "MMFFGetMoleculeProperties",
    },

    # ── Reactions ──
    "reaction": {
        "AllChem.ReactionFromSmarts", "ReactionFromSmarts",
        "rdChemReactions.ReactionFromSmarts",
        "AllChem.ReactionFromRxnFile", "ReactionFromRxnFile",
        "RunReactants", "rxn.RunReactants",
    },

    # ── Scaffolds ──
    "scaffold": {
        "MurckoScaffold.GetScaffoldForMol", "GetScaffoldForMol",
        "MurckoScaffold.MakeScaffoldGeneric", "MakeScaffoldGeneric",
    },

    # ── Salt/Fragment removal ──
    "salt_strip": {
        "SaltRemover", "rdmolops.GetMolFrags",
        "Chem.GetMolFrags", "GetMolFrags",
    },

    # ── Standardization ──
    "standardize": {
        "MolStandardize", "rdMolStandardize",
        "Uncharger", "LargestFragmentChooser",
        "Normalizer", "TautomerCanonicalizer",
    },

    # ── MCS ──
    "mcs": {
        "rdFMCS.FindMCS", "FindMCS",
    },

    # ── R-group decomposition ──
    "rgroup": {
        "rdRGroupDecomposition.RGroupDecompose", "RGroupDecompose",
        "rdRGroupDecomposition.RGroupDecomposition", "RGroupDecomposition",
    },

    # ── BRICS ──
    "brics": {
        "BRICS.BRICSDecompose", "BRICSDecompose",
        "BRICS.BRICSBuild", "BRICSBuild",
    },

    # ── Atom/Bond properties ──
    "atom_props": {
        "GetAtomicNum", "GetSymbol", "GetFormalCharge",
        "GetNumImplicitHs", "GetNumExplicitHs",
        "GetIsAromatic", "GetHybridization", "GetDegree",
        "GetTotalNumHs", "GetChiralTag",
    },
    "bond_props": {
        "GetBondType", "GetBondTypeAsDouble", "GetIsAromatic",
        "GetIsConjugated", "GetStereo",
    },

    # ── Ring info ──
    "ring_info": {
        "GetRingInfo", "GetSymmSSSR", "GetSSSR",
        "mol.GetRingInfo",
    },

    # ── Aromaticity / Kekulize ──
    "aromaticity": {
        "Chem.Kekulize", "Kekulize",
        "Chem.SetAromaticity", "SetAromaticity",
    },

    # ── AddHs / RemoveHs ──
    "hydrogens": {
        "Chem.AddHs", "AddHs",
        "Chem.RemoveHs", "RemoveHs",
    },

    # ── Compute2DCoords ──
    "coords_2d": {
        "AllChem.Compute2DCoords", "Compute2DCoords",
        "rdDepictor.Compute2DCoords",
    },

    # ── Crippen ──
    "crippen": {
        "Crippen.MolLogP", "Crippen.MolMR",
        "rdMolDescriptors.CalcCrippenDescriptors", "CalcCrippenDescriptors",
    },

    # ── Pandas / ML ──
    "pandastools": {
        "PandasTools", "PandasTools.LoadSDF", "PandasTools.WriteSDF",
        "PandasTools.AddMoleculeColumnToFrame",
    },
}


# ============================================================
# AST-based API call extractor
# ============================================================

class _APICallVisitor(ast.NodeVisitor):
    """Walk AST and collect all attribute-style calls like Chem.MolFromSmiles."""

    def __init__(self):
        self.calls: Set[str] = set()
        self.imports: Dict[str, str] = {}  # alias -> module

    def visit_ImportFrom(self, node: ast.ImportFrom):
        """Track 'from X import Y as Z' to resolve aliases."""
        if node.module and node.names:
            for alias in node.names:
                name = alias.asname or alias.name
                self.imports[name] = f"{node.module}.{alias.name}"
        self.generic_visit(node)

    def visit_Import(self, node: ast.Import):
        for alias in node.names:
            name = alias.asname or alias.name
            self.imports[name] = alias.name
        self.generic_visit(node)

    def visit_Call(self, node: ast.Call):
        """Extract the qualified name of every function/method call."""
        name = self._resolve_call_name(node.func)
        if name:
            self.calls.add(name)
        self.generic_visit(node)

    def visit_Attribute(self, node: ast.Attribute):
        """Also capture method-style access even without call (e.g. mol.GetRingInfo)."""
        name = self._get_attr_chain(node)
        if name:
            self.calls.add(name)
        self.generic_visit(node)

    def _resolve_call_name(self, node) -> Optional[str]:
        if isinstance(node, ast.Name):
            return node.id
        elif isinstance(node, ast.Attribute):
            return self._get_attr_chain(node)
        return None

    def _get_attr_chain(self, node, max_depth=4) -> Optional[str]:
        """Reconstruct dotted name like 'AllChem.GetMorganFingerprintAsBitVect'."""
        parts = []
        current = node
        depth = 0
        while isinstance(current, ast.Attribute) and depth < max_depth:
            parts.append(current.attr)
            current = current.value
            depth += 1
        if isinstance(current, ast.Name):
            parts.append(current.id)
        parts.reverse()
        return ".".join(parts) if parts else None


def extract_api_calls(source_code: str) -> Set[str]:
    """
    Parse Python source and return all function/method call names found.
    Returns a set of dotted names like {'Chem.MolFromSmiles', 'AllChem.EmbedMolecule'}.
    """
    try:
        tree = ast.parse(source_code)
    except SyntaxError:
        return set()

    visitor = _APICallVisitor()
    visitor.visit(tree)
    return visitor.calls


def extract_api_calls_from_file(filepath: str) -> Set[str]:
    """Extract API calls from a Python file."""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            source = f.read()
        return extract_api_calls(source)
    except Exception:
        return set()


# ============================================================
# Categorize API calls
# ============================================================

def categorize_calls(calls: Set[str]) -> Set[str]:
    """
    Map raw API call names to functional categories.
    Returns a set of category names like {'mol_from_smiles', 'fp_morgan', 'similarity'}.
    """
    categories = set()
    for call in calls:
        # Check each part of the dotted name against category entries
        # e.g. "AllChem.GetMorganFingerprintAsBitVect" should match
        # even if stored as just "GetMorganFingerprintAsBitVect"
        for cat_name, cat_apis in RDKIT_API_CATEGORIES.items():
            for api in cat_apis:
                # Match if the call ends with the API name, or contains it
                if call == api or call.endswith("." + api) or api.endswith("." + call.split(".")[-1]):
                    categories.add(cat_name)
                    break
                # Also match the last segment (method name)
                call_method = call.split(".")[-1]
                api_method = api.split(".")[-1]
                if call_method == api_method and call_method not in (
                    # Exclude overly generic names
                    "append", "keys", "values", "items", "get", "set",
                    "read", "write", "open", "close", "print", "len",
                    "range", "enumerate", "zip", "map", "filter", "sorted",
                    "round", "abs", "max", "min", "sum", "any", "all",
                    "isinstance", "type", "str", "int", "float", "bool", "list", "dict", "tuple", "set",
                ):
                    categories.add(cat_name)
                    break
    return categories


# ============================================================
# Key API extraction (solution-focused)
# ============================================================

# These categories are "infrastructure" — almost every solution uses them.
# They should NOT count as "key" differentiating APIs.
INFRASTRUCTURE_CATEGORIES = {
    "mol_from_smiles",
    "mol_from_smarts",
    "hydrogens",  # AddHs/RemoveHs is often a preprocessing step
}

# Functional equivalence groups: categories within the same group
# are considered interchangeable (they achieve the same goal via different APIs).
# When comparing, if sol uses category A and pred uses category B,
# and A and B are in the same equivalence group, it counts as a match.
EQUIVALENT_GROUPS = [
    # All visualization approaches are equivalent
    {"mol_draw_image", "mol_draw_svg", "draw_general"},
    # Different fingerprint types are functionally equivalent
    {"fp_morgan", "fp_rdkit", "fp_maccs", "fp_topological"},
    # Different force fields for optimization
    {"conformer_optimize", "conformer_energy"},
    # MolBlock and PDB output are both 3D structure representations
    {"mol_to_molblock", "mol_to_pdb"},
    # Different I/O formats for reading molecules
    {"mol_from_molblock", "mol_from_pdb", "mol_from_mol2", "sdf_io"},
    # SMILES and SMARTS output
    {"mol_to_smiles", "mol_to_smarts"},
    # Crippen LogP vs Descriptors LogP
    {"crippen", "descriptors_basic"},
]

def _build_equivalence_map() -> Dict[str, int]:
    """Build a mapping from category -> group_id for fast lookup."""
    eq_map = {}
    for group_id, group in enumerate(EQUIVALENT_GROUPS):
        for cat in group:
            eq_map[cat] = group_id
    return eq_map

_EQUIVALENCE_MAP = _build_equivalence_map()


def _normalize_categories(categories: Set[str]) -> Set[str]:
    """
    Normalize categories by replacing each with its equivalence group representative.
    This way, 'mol_draw_image' and 'mol_draw_svg' both become the same normalized key.
    """
    normalized = set()
    for cat in categories:
        if cat in _EQUIVALENCE_MAP:
            # Use the group_id as a canonical representative
            group_id = _EQUIVALENCE_MAP[cat]
            # Use the first element of the group as canonical name
            canonical = sorted(EQUIVALENT_GROUPS[group_id])[0]
            normalized.add(canonical)
        else:
            normalized.add(cat)
    return normalized


def extract_key_categories(filepath: str) -> Set[str]:
    """
    Extract the key (non-infrastructure) API categories from a solution file.
    These represent the "core work" the solution does.
    """
    calls = extract_api_calls_from_file(filepath)
    categories = categorize_calls(calls)
    # Remove infrastructure categories
    key_cats = categories - INFRASTRUCTURE_CATEGORIES
    return key_cats


# ============================================================
# API match comparison
# ============================================================

def compare_api_match(
    solution_path: str,
    prediction_path: str,
    min_overlap_ratio: float = 0.5,
) -> Tuple[bool, str]:
    """
    Check if prediction code calls a sufficient subset of the solution's
    key RDKit APIs, accounting for functional equivalence.

    Args:
        solution_path: Path to solution .py file
        prediction_path: Path to prediction .py file
        min_overlap_ratio: Minimum fraction of solution's key categories
                          that must appear in prediction (default 0.5 = 50%)

    Returns:
        (is_match, detail_string)
    """
    sol_cats = extract_key_categories(solution_path)
    pred_cats = extract_key_categories(prediction_path)

    if not sol_cats:
        # Solution has no distinguishing API calls — can't evaluate this way
        return False, "no key APIs in solution"

    # Normalize both sides using equivalence groups
    sol_norm = _normalize_categories(sol_cats)
    pred_norm = _normalize_categories(pred_cats)

    overlap = sol_norm & pred_norm
    overlap_ratio = len(overlap) / len(sol_norm)

    detail = (
        f"sol_apis={sorted(sol_cats)}, "
        f"pred_apis={sorted(pred_cats)}, "
        f"sol_norm={sorted(sol_norm)}, "
        f"pred_norm={sorted(pred_norm)}, "
        f"overlap={sorted(overlap)}, "
        f"ratio={overlap_ratio:.2f}"
    )

    if overlap_ratio >= min_overlap_ratio:
        return True, detail
    else:
        return False, detail


# ============================================================
# Combined fallback evaluation
# ============================================================

def api_fallback_check(
    solution_path: str,
    prediction_path: str,
    pred_executed_ok: bool,
    min_overlap_ratio: float = 0.5,
) -> Tuple[bool, str]:
    """
    Fallback evaluation when output comparison fails.

    Conditions for a pass:
        1. Prediction code executed without error (pred_executed_ok=True)
        2. Prediction code calls >= min_overlap_ratio of solution's key APIs

    Args:
        solution_path: Path to solution .py file
        prediction_path: Path to prediction .py file
        pred_executed_ok: Whether the prediction executed without errors
        min_overlap_ratio: Minimum API overlap ratio

    Returns:
        (is_pass, detail_string)
    """
    if not pred_executed_ok:
        return False, "prediction did not execute successfully"

    if not os.path.exists(prediction_path):
        return False, "prediction file not found"

    api_match, api_detail = compare_api_match(
        solution_path, prediction_path, min_overlap_ratio
    )

    if api_match:
        return True, f"api_fallback_pass: {api_detail}"
    else:
        return False, f"api_fallback_fail: {api_detail}"
