from rdkit import Chem
import random


def level_function(mol):
    """给定分子，随机删除一个原子。"""
    try:
        mol = Chem.MolFromSmiles(mol)
        rwmol = Chem.RWMol(mol)
        num_atoms = rwmol.GetNumAtoms()
        if num_atoms <= 1:
            return None
        # Find removable atoms (not in rings to avoid breaking ring structure)
        removable = []
        for idx in range(num_atoms):
            atom = rwmol.GetAtomWithIdx(idx)
            if not atom.IsInRing() and atom.GetDegree() <= 1:
                removable.append(idx)
        if not removable:
            # If no safe removable atoms, pick any non-ring atom
            for idx in range(num_atoms):
                atom = rwmol.GetAtomWithIdx(idx)
                if not atom.IsInRing():
                    removable.append(idx)
        if not removable:
            removable = list(range(num_atoms))
        idx_to_remove = random.choice(removable)
        rwmol.RemoveAtom(idx_to_remove)
        try:
            Chem.SanitizeMol(rwmol)
            return Chem.MolToSmiles(rwmol)
        except Exception:
            return Chem.MolToSmiles(rwmol)
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    result = level_function("c1ccccc1CCO")
    print(f"删除一个原子后的分子: {result}")
