from rdkit import Chem
import random


def level_function(mol):
    """给定分子，随机交换两个原子。"""
    try:
        mol = Chem.MolFromSmiles(mol)
        rwmol = Chem.RWMol(mol)
        num_atoms = rwmol.GetNumAtoms()
        if num_atoms < 2:
            return None
        indices = list(range(num_atoms))
        idx1, idx2 = random.sample(indices, 2)
        atom1 = rwmol.GetAtomWithIdx(idx1)
        atom2 = rwmol.GetAtomWithIdx(idx2)
        # Swap atomic numbers
        anum1 = atom1.GetAtomicNum()
        anum2 = atom2.GetAtomicNum()
        atom1.SetAtomicNum(anum2)
        atom2.SetAtomicNum(anum1)
        try:
            Chem.SanitizeMol(rwmol)
            return Chem.MolToSmiles(rwmol)
        except Exception:
            # Revert if sanitization fails
            atom1.SetAtomicNum(anum1)
            atom2.SetAtomicNum(anum2)
            return Chem.MolToSmiles(mol)
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    result = level_function("c1ccccc1CCO")
    print(f"交换两个原子后的分子: {result}")
