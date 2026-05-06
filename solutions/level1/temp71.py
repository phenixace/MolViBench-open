from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        return rdMolDescriptors.CalcNumAromaticRings(mol_obj)
    except Exception as e:
        print(e)
        return None
