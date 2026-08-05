from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol):



    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        rxn = AllChem.ReactionFromSmarts('[c:1][OH]>>[c:1]OC')
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
    smiles = 'Oc1ccccc1'
    result = level_function(smiles)
    print(f'Output: {result}')
