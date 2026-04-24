from rdkit import Chem
import random


def level_function(mol):
    """随机扰动分子结构（原子级别）。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        rw = Chem.RWMol(mol_obj)
        num_atoms = rw.GetNumAtoms()
        if num_atoms == 0:
            return None

        # 随机选择一个非氢原子，将其原子序数替换为另一个常见元素
        replacements = [6, 7, 8, 9, 16, 17]  # C, N, O, F, S, Cl
        idx = random.randint(0, num_atoms - 1)
        orig_num = rw.GetAtomWithIdx(idx).GetAtomicNum()
        candidates = [z for z in replacements if z != orig_num]
        if not candidates:
            return None

        new_num = random.choice(candidates)
        rw.GetAtomWithIdx(idx).SetAtomicNum(new_num)

        try:
            Chem.SanitizeMol(rw)
            return Chem.MolToSmiles(rw)
        except Exception:
            return None
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccccc1"
    for _ in range(5):
        print(f"原子级扰动: {level_function(smiles)}")
