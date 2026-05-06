from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        fp = AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(mol_obj, nBits=2048)
        return fp.ToBitString()
    except Exception as e:
        print(e)
        return None
