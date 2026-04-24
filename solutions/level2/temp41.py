from rdkit import Chem
import random


def level_function(mol):
    """给定分子，随机添加一个甲基。"""
    try:
        mol = Chem.MolFromSmiles(mol)
        rwmol = Chem.RWMol(mol)
        # Find atoms with available valence
        candidates = []
        for atom in rwmol.GetAtoms():
            default_valence = Chem.GetPeriodicTable().GetDefaultValence(atom.GetAtomicNum())
            if isinstance(default_valence, tuple):
                max_valence = max(default_valence)
            else:
                max_valence = default_valence
            current_valence = atom.GetTotalValence()
            if current_valence < max_valence:
                candidates.append(atom.GetIdx())
        if not candidates:
            return Chem.MolToSmiles(mol)
        target_idx = random.choice(candidates)
        new_idx = rwmol.AddAtom(Chem.Atom(6))  # Carbon
        rwmol.AddBond(target_idx, new_idx, Chem.rdchem.BondType.SINGLE)
        Chem.SanitizeMol(rwmol)
        return Chem.MolToSmiles(rwmol)
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    result = level_function("c1ccccc1")
    print(f"添加甲基后的分子: {result}")
