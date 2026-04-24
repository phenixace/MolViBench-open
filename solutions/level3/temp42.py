from rdkit import Chem
from rdkit.Chem import AllChem, BRICS
import random


def level_function(mol):
    """随机扰动分子结构（侧链级别）。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # 使用 BRICS 分解来识别侧链切割点
        fragments = list(BRICS.BRICSDecompose(mol_obj))
        if len(fragments) < 2:
            # 如果无法 BRICS 分解, 尝试随机删除一个末端原子(侧链)
            rw = Chem.RWMol(mol_obj)
            terminal_atoms = [
                atom.GetIdx() for atom in rw.GetAtoms()
                if atom.GetDegree() == 1 and atom.GetAtomicNum() != 1
            ]
            if not terminal_atoms:
                return None
            remove_idx = random.choice(terminal_atoms)
            rw.RemoveAtom(remove_idx)
            try:
                Chem.SanitizeMol(rw)
                return Chem.MolToSmiles(rw)
            except Exception:
                return None

        # 随机替换一个片段
        substituents = ["C", "CC", "O", "N", "F", "Cl"]
        sub = random.choice(substituents)
        sub_mol = Chem.MolFromSmiles(sub)
        if sub_mol is None:
            return None

        # 用反应方式替换侧链: 删除末端原子并添加新原子
        rw = Chem.RWMol(mol_obj)
        terminal_atoms = [
            atom.GetIdx() for atom in rw.GetAtoms()
            if atom.GetDegree() == 1 and atom.GetAtomicNum() != 1
        ]
        if not terminal_atoms:
            return None

        remove_idx = random.choice(terminal_atoms)
        neighbor = rw.GetAtomWithIdx(remove_idx).GetNeighbors()[0].GetIdx()
        rw.RemoveAtom(remove_idx)
        # 添加新的取代基
        new_idx = rw.AddAtom(Chem.Atom(sub_mol.GetAtomWithIdx(0).GetAtomicNum()))
        # 调整 neighbor index (如果 remove_idx < neighbor)
        adj_neighbor = neighbor if remove_idx > neighbor else neighbor
        rw.AddBond(adj_neighbor, new_idx, Chem.BondType.SINGLE)

        try:
            Chem.SanitizeMol(rw)
            return Chem.MolToSmiles(rw)
        except Exception:
            return None
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"  # 布洛芬
    for _ in range(5):
        print(f"侧链级扰动: {level_function(smiles)}")
