from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

def level_function(mol1, mol2):
    try:
        mol1 = Chem.MolFromSmiles(mol1)
        if mol1 is None:
            return None
        mol2 = Chem.MolFromSmiles(mol2)
        if mol2 is None:
            return None
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
        return DataStructs.DiceSimilarity(fp1, fp2)
    except Exception as e:
        print(e)
        return None
