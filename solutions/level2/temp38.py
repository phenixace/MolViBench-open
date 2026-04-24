from rdkit import Chem
import random


def level_function(mol):
    """给定分子，随机删除一个支链。"""
    try:
        mol = Chem.MolFromSmiles(mol)
        rwmol = Chem.RWMol(mol)
        # Find terminal atoms (degree 1, not in ring) as side chain endpoints
        terminal_atoms = []
        for atom in rwmol.GetAtoms():
            if atom.GetDegree() == 1 and not atom.IsInRing():
                terminal_atoms.append(atom.GetIdx())
        if not terminal_atoms:
            return Chem.MolToSmiles(mol)
        # Pick a random terminal atom and trace back the side chain
        start_idx = random.choice(terminal_atoms)
        to_remove = []
        current = rwmol.GetAtomWithIdx(start_idx)
        while current.GetDegree() <= 1 and not current.IsInRing():
            to_remove.append(current.GetIdx())
            neighbors = [n for n in current.GetNeighbors() if n.GetIdx() not in to_remove]
            if not neighbors:
                break
            current = neighbors[0]
            if current.GetDegree() > 2 or current.IsInRing():
                break
            # Continue tracing if atom has degree <= 2 and is not in ring
        # Remove atoms in reverse order of index to maintain indices
        for idx in sorted(to_remove, reverse=True):
            rwmol.RemoveAtom(idx)
        Chem.SanitizeMol(rwmol)
        return Chem.MolToSmiles(rwmol)
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    result = level_function("c1ccccc1CCO")
    print(f"删除支链后的分子: {result}")
