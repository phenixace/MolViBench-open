from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol):
    """
    生成分子的 Morgan 指纹（半径=3）。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=2048)
        return list(fp)
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CCO"
    fp = level_function(smiles)
    print(f"Morgan 指纹（半径=3）长度: {len(fp)}, 前20位: {fp[:20]}")
