from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(inhibitor_smiles, cysteine_smiles="NC(CS)C(=O)O"):
    """模拟共价抑制剂（含丙烯酰胺弹头）与半胱氨酸的 Michael 加成。"""
    try:
        inhibitor = Chem.MolFromSmiles(inhibitor_smiles)
        cysteine = Chem.MolFromSmiles(cysteine_smiles)
        if inhibitor is None or cysteine is None:
            return None

        # Check if inhibitor contains acrylamide warhead
        acrylamide_pattern = Chem.MolFromSmarts("[#6]C(=O)NC=C")
        if not inhibitor.HasSubstructMatch(acrylamide_pattern):
            # Also check simple vinyl amide
            acrylamide_pattern2 = Chem.MolFromSmarts("C=CC(=O)N")
            if not inhibitor.HasSubstructMatch(acrylamide_pattern2):
                return {"has_warhead": False, "product": None}

        # Michael addition: Cys-SH + C=C-C(=O)-N -> Cys-S-CH2-CH2-C(=O)-N
        rxn = AllChem.ReactionFromSmarts(
            "[S:1][H].[CH2:2]=[CH:3][C:4](=[O:5])[N:6]>>[S:1][CH2:2][CH2:3][C:4](=[O:5])[N:6]"
        )

        products = rxn.RunReactants((cysteine, inhibitor))
        if not products:
            return {"has_warhead": True, "product": None}

        results = set()
        for prod_set in products:
            for prod in prod_set:
                try:
                    Chem.SanitizeMol(prod)
                    results.add(Chem.MolToSmiles(prod))
                except Exception:
                    pass

        return {
            "has_warhead": True,
            "product": sorted(list(results)) if results else None
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    # Acrylamide warhead: simple example
    smiles = "C=CC(=O)Nc1ccccc1"
    result = level_function(smiles)
    if result:
        print(f"Has warhead: {result['has_warhead']}")
        print(f"Michael addition product: {result['product']}")
