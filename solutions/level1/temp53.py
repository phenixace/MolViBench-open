from rdkit import Chem
from rdkit.Chem import Descriptors


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        Chem.AssignStereochemistry(mol_obj, cleanIt=True, force=True)
        chiral_centers = Chem.FindMolChiralCenters(mol_obj, includeUnassigned=True)

        result = [(idx, label) for idx, label in chiral_centers]
        return result
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'C[C@H](O)[C@@H](F)Cl'
    print(f'Output: {level_function(smiles)}')
