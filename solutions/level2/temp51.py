from rdkit import Chem
from rdkit.Chem import MACCSkeys


def level_function(mol):
    """生成分子的 MACCS Keys 指纹。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        fp = MACCSkeys.GenMACCSKeys(mol_obj)
        return fp.ToBitString()
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccccc1"
    result = level_function(smiles)
    print(f"MACCS Keys 指纹 (长度={len(result) if result else 0}): {result}")
