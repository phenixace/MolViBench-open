from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

def level_function(smiles1: str, smiles2: str,
                         radius: int = 2,
                         nBits: int = 2048,
                         useFeatures: bool = False) -> float:
    """
    计算两条 SMILES 之间的 Tanimoto 相似度。

    参数说明
    ----------
    smiles1, smiles2 : str
        输入的 SMILES 字符串。
    radius : int
        Morgan 指纹的半径，默认 2（对应 6‑隧道）。常见取值 2~3。
    nBits : int
        指纹长度，默认 2048。
    useFeatures : bool
        是否使用 RDKit 的化学特征指纹（RDKFingerprint）。若 False 使用 Morgan 指纹。

    返回
    ------
    similarity : float
        0~1 之间的 Tanimoto 相似度。
    """
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

# 示例
if __name__ == "__main__":
    smiles_a = "CC1=CC=CC=C1"        # 甲苯
    smiles_b = "CC1=CC=CC=C1C"      # 乙基苯

    sim = tanimoto_from_smiles(smiles_a, smiles_b)
    print(f"Tanimoto 相似度 (Morgan, radius=2, 2048 bits): {sim:.4f}")
