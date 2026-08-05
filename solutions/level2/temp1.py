from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

def level_function(smiles1: str, smiles2: str,
                         radius: int = 2,
                         nBits: int = 2048,
                         useFeatures: bool = False) -> float:



















    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    if mol1 is None or mol2 is None:
        raise ValueError("One of the SMILES could not be parsed.")

    if useFeatures:
        fp1 = Chem.RDKFingerprint(mol1, fpSize=nBits)
        fp2 = Chem.RDKFingerprint(mol2, fpSize=nBits)
    else:
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius, nBits=nBits)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius, nBits=nBits)

    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    return similarity


if __name__ == '__main__':
    smiles_a = 'CC1=CC=CC=C1'
    smiles_b = 'CC1=CC=CC=C1C'
    sim = level_function(smiles_a, smiles_b)
    print(f'Output: {sim:.4f}')
