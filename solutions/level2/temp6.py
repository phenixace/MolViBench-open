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
        distance = 0
        for i in range(fp1.GetNumBits()):
            if fp1.GetBit(i) != fp2.GetBit(i):
                distance += 1
        return distance
    except Exception as e:
        print(e)
        return None
