from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol):
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=2048)
        return list(fp)
    except Exception as e:
        print(e)
        return None
