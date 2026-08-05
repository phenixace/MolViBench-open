import selfies as sf
from rdkit import Chem


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        canonical = Chem.MolToSmiles(mol_obj)
        selfies_str = sf.encoder(canonical)
        return selfies_str
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'CCO'
    print(f'Output: {level_function(smiles)}')
