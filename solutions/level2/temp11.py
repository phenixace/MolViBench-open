from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol, substituent="C"):



    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        rxn_smarts = f'[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1>>[c:1]({substituent})1[cH:2][c:3]({substituent})[cH:4][cH:5][cH:6]1'
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        products = rxn.RunReactants((mol_obj,))
        unique_smiles = set()
        for product_tuple in products:
            for product in product_tuple:
                Chem.SanitizeMol(product)
                smi = Chem.MolToSmiles(product)
                unique_smiles.add(smi)
        if unique_smiles:
            return sorted(list(unique_smiles))[0]
        return None
    except Exception as e:
        print(e)
        return None

if __name__ == '__main__':
    smiles = 'c1ccccc1'
    result = level_function(smiles, 'C')
    print(f'Output: {result}')
