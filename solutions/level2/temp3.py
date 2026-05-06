from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

def level_function(mol1, mol2):
    try:
        m1 = Chem.MolFromSmiles(mol1)
        if m1 is None:
            return None
        m2 = Chem.MolFromSmiles(mol2)
        if m2 is None:
            return None
        fp1 = AllChem.GetMorganFingerprintAsBitVect(m1, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(m2, 2, nBits=2048)
        return DataStructs.CosineSimilarity(fp1, fp2)
    except Exception as e:
        print(e)
        return None
