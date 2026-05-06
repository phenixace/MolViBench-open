from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol1, mol2):
    try:
        reaction_smarts = '[C:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2]=[C:3][C:4][C:6][C:5]1'
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        diene = Chem.MolFromSmiles(mol1)
        dienophile = Chem.MolFromSmiles(mol2)
        if diene is None or dienophile is None:
            return None
        products = rxn.RunReactants((diene, dienophile))
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
