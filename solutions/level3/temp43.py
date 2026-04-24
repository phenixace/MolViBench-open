from rdkit import Chem
import random


def level_function(mol):
    """随机扰动分子结构（环级别）。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        ring_info = mol_obj.GetRingInfo()
        atom_rings = ring_info.AtomRings()
        if not atom_rings:
            return None

        # 随机选择一个环
        ring = list(random.choice(list(atom_rings)))

        rw = Chem.RWMol(mol_obj)

        # 策略: 随机修改环中的一个原子
        replacements = [6, 7, 8, 16]  # C, N, O, S
        target_idx = random.choice(ring)
        orig_num = rw.GetAtomWithIdx(target_idx).GetAtomicNum()
        candidates = [z for z in replacements if z != orig_num]
        if not candidates:
            return None

        new_num = random.choice(candidates)
        rw.GetAtomWithIdx(target_idx).SetAtomicNum(new_num)

        try:
            Chem.SanitizeMol(rw)
            return Chem.MolToSmiles(rw)
        except Exception:
            return None
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccccc1"  # 苯
    for _ in range(5):
        print(f"环级扰动: {level_function(smiles)}")
