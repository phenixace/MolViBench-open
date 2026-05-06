from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol1, mol2):
    try:
        reaction_smarts = '[C:1](=[O:2])[OH:3].[NH2:4][C:5]>>[C:1](=[O:2])[NH:4][C:5]'
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        acid = Chem.MolFromSmiles(mol1)
        amine = Chem.MolFromSmiles(mol2)
        if acid is None or amine is None:
            return None
        products = rxn.RunReactants((acid, amine))
        if not products:
            return None
        result_smiles = []
        for product_set in products:
            for product in product_set:
                Chem.SanitizeMol(product)
                smi = Chem.MolToSmiles(product)
                if smi not in result_smiles:
                    result_smiles.append(smi)
        return result_smiles
    except Exception as e:
        print(e)
        return None
