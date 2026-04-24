from rdkit import Chem
from rdkit.Chem import AllChem
import random

def level_function(mol):
    """
    在分子中随机替换一个氢为碘。
    """
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        # 添加显式氢
        mol_h = Chem.AddHs(mol_obj)
        # 找到所有氢原子的索引
        h_indices = [atom.GetIdx() for atom in mol_h.GetAtoms() if atom.GetAtomicNum() == 1]
        if not h_indices:
            return None
        # 随机选择一个氢原子
        h_idx = random.choice(h_indices)
        # 使用 RWMol 替换氢为碘
        rw_mol = Chem.RWMol(mol_h)
        rw_mol.GetAtomWithIdx(h_idx).SetAtomicNum(53)  # 53 = I
        Chem.SanitizeMol(rw_mol)
        rw_mol = Chem.RemoveHs(rw_mol)
        return Chem.MolToSmiles(rw_mol)
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "c1ccccc1"  # 苯
    result = level_function(smiles)
    print(f"随机替换一个氢为碘: {result}")
