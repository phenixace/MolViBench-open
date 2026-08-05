from rdkit import Chem
from rdkit.Chem import Descriptors


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        heavy = mol_obj.GetNumHeavyAtoms()
        if heavy == 0:
            return 0.0
        aromatic = sum(1 for atom in mol_obj.GetAtoms() if atom.GetIsAromatic())
        return round(aromatic / heavy, 4)
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccc(CC)cc1'
    print(f'Output: {level_function(smiles)}')
