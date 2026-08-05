from rdkit import Chem
from rdkit.Chem import GraphDescriptors


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        j = GraphDescriptors.BalabanJ(mol_obj)
        return round(j, 4)
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccccc1'
    print(f'Output: {level_function(smiles)}')
