from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(azide_smiles, alkyne_smiles):
    try:
        azide = Chem.MolFromSmiles(azide_smiles)
        alkyne = Chem.MolFromSmiles(alkyne_smiles)
        if azide is None or alkyne is None:
            return None

        rxn_smarts = "[C:1]#[CH1:2].[N:3]=[N:4]=[N:5]>>[C:1]1=NN([N:3])[N:5]=[C:2]1"
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)

        products = rxn.RunReactants((alkyne, azide))
        if not products:
            products = rxn.RunReactants((azide, alkyne))

        if not products:
            return None

        results = set()
        for prod_set in products:
            for prod in prod_set:
                try:
                    Chem.SanitizeMol(prod)
                    results.add(Chem.MolToSmiles(prod))
                except Exception:
                    pass

        return sorted(list(results)) if results else None
    except Exception as e:
        print(e)
        return None
