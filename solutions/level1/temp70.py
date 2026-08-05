from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        return rdMolDescriptors.CalcNumAliphaticRings(mol_obj)
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'C1CCC(CC1)c1ccccc1'
    print(f'Output: {level_function(smiles)}')
