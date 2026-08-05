from rdkit import Chem
import random


def level_function(mol):

    try:
        mol = Chem.MolFromSmiles(mol)
        rwmol = Chem.RWMol(mol)

        terminal_atoms = []
        for atom in rwmol.GetAtoms():
            if atom.GetDegree() == 1 and not atom.IsInRing():
                terminal_atoms.append(atom.GetIdx())
        if not terminal_atoms:
            return Chem.MolToSmiles(mol)

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


        for idx in sorted(to_remove, reverse=True):
            rwmol.RemoveAtom(idx)
        Chem.SanitizeMol(rwmol)
        return Chem.MolToSmiles(rwmol)
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    result = level_function('c1ccccc1CCO')
    print(f'Output: {result}')
