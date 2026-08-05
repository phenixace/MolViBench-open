from rdkit import Chem
from rdkit.Chem import Descriptors


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        return Descriptors.NHOHCount(mol_obj)
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'CC(O)c1ccc(N)cc1'
    print(f'Output: {level_function(smiles)}')
