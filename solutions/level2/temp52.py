from rdkit import Chem
from rdkit.Chem import RDKFingerprint

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        fp = RDKFingerprint(mol_obj)
        return fp.ToBitString()
    except Exception as e:
        print(e)
        return None
