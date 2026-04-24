from rdkit import Chem
from rdkit.Chem import RDKFingerprint


def level_function(mol):
    """生成分子的 RDKit 拓扑指纹。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        fp = RDKFingerprint(mol_obj)
        return fp.ToBitString()
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccccc1"
    result = level_function(smiles)
    print(f"RDKit 拓扑指纹 (长度={len(result) if result else 0})")
