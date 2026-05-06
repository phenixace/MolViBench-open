from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(inhibitor_smiles, cysteine_smiles="NC(CS)C(=O)O"):
    try:
        inhibitor = Chem.MolFromSmiles(inhibitor_smiles)
        cysteine = Chem.MolFromSmiles(cysteine_smiles)
        if inhibitor is None or cysteine is None:
            return None

        acrylamide_pattern = Chem.MolFromSmarts("[#6]C(=O)NC=C")
        if not inhibitor.HasSubstructMatch(acrylamide_pattern):
            acrylamide_pattern2 = Chem.MolFromSmarts("C=CC(=O)N")
            if not inhibitor.HasSubstructMatch(acrylamide_pattern2):
                return {"has_warhead": False, "product": None}

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
