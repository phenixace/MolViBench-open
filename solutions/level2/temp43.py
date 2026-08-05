from rdkit import Chem
import random


def level_function(mol):

    try:
        mol = Chem.MolFromSmiles(mol)
        rwmol = Chem.RWMol(mol)

        candidates = []
        for atom in rwmol.GetAtoms():
            default_valence = Chem.GetPeriodicTable().GetDefaultValence(atom.GetAtomicNum())
            if isinstance(default_valence, tuple):
                max_valence = max(default_valence)
            else:
                max_valence = default_valence
            current_valence = atom.GetTotalValence() - atom.GetNumImplicitHs()
            if current_valence < max_valence:
                candidates.append(atom.GetIdx())
        if not candidates:
            return Chem.MolToSmiles(mol)
        target_idx = random.choice(candidates)

        halogen = random.choice([9, 17, 35, 53])
        new_idx = rwmol.AddAtom(Chem.Atom(halogen))
        rwmol.AddBond(target_idx, new_idx, Chem.rdchem.BondType.SINGLE)
        Chem.SanitizeMol(rwmol)
        return Chem.MolToSmiles(rwmol)
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    result = level_function('c1ccccc1')
    print(f'Output: {result}')
