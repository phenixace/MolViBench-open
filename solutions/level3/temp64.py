from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(inhibitor_smiles, cysteine_smiles="NC(CS)C(=O)O"):
    try:
        inhibitor = Chem.MolFromSmiles(inhibitor_smiles)
        cysteine = Chem.MolFromSmiles(cysteine_smiles)
        if inhibitor is None or cysteine is None:
            return None



        acrylamide_pattern = Chem.MolFromSmarts("C=CC(=O)N")
        if not inhibitor.HasSubstructMatch(acrylamide_pattern):
            return {"has_warhead": False, "product": None}





        rxn_smarts = "[S:1].[CH2:2]=[CH:3][C:4]=[O:5]>>[S:1][CH2:2][CH2:3][C:4]=[O:5]"
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)


        products = rxn.RunReactants((cysteine, inhibitor))

        if not products:
            return {"has_warhead": True, "product": None}

        results = set()
        for prod_set in products:
            for prod in prod_set:
                try:

                    prod.UpdatePropertyCache()
                    Chem.SanitizeMol(prod)
                    results.add(Chem.MolToSmiles(prod))
                except Exception:
                    continue

        return {
            "has_warhead": True,
            "product": sorted(list(results)) if results else None
        }
    except Exception as e:
        print(f"Error: {e}")
        return None

if __name__ == '__main__':
    smiles = 'C=CC(=O)Nc1ccccc1'
    result = level_function(smiles)
    if result:
        print(f"Output: {result['has_warhead']}")
        print(f"Output: {result['product']}")
