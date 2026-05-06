from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol1, mol2):
    try:
        acid = Chem.MolFromSmiles(mol1)
        amine = Chem.MolFromSmiles(mol2)
        if acid is None or amine is None:
            return None

        rxn_smarts = '[C:1](=O)[OH].[N:2]>>[C:1](=O)[N:2]'
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        products = rxn.RunReactants((acid, amine))

        product_smiles = []
        for prod_set in products:
            for prod in prod_set:
                Chem.SanitizeMol(prod)
                smi = Chem.MolToSmiles(prod)
                if smi not in product_smiles:
                    product_smiles.append(smi)
        return product_smiles
    except Exception as e:
        print(e)
        return None
