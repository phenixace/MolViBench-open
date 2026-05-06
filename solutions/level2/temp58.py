from rdkit import Chem
from rdkit.Chem import BRICS

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        fragments = BRICS.BRICSDecompose(mol_obj)
        return sorted(list(fragments))
    except Exception as e:
        print(e)
        return None
