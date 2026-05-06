from rdkit import Chem
from rdkit.Chem import Descriptors

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        exact_mw = Descriptors.ExactMolWt(mol_obj)
        return round(exact_mw, 4)
    except Exception as e:
        print(e)
        return None
