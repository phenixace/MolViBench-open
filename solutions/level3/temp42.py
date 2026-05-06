from rdkit import Chem
from rdkit.Chem import AllChem, BRICS
import random

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        fragments = list(BRICS.BRICSDecompose(mol_obj))
        if len(fragments) < 2:
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

        substituents = ["C", "CC", "O", "N", "F", "Cl"]
        sub = random.choice(substituents)
        sub_mol = Chem.MolFromSmiles(sub)
        if sub_mol is None:
            return None

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
        new_idx = rw.AddAtom(Chem.Atom(sub_mol.GetAtomWithIdx(0).GetAtomicNum()))
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
