from rdkit import Chem
from rdkit.Chem import Descriptors


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        fsp3 = Descriptors.FractionCSP3(mol_obj)
        return round(fsp3, 4)
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'CCCC(=O)O'
    print(f'Output: {level_function(smiles)}')
