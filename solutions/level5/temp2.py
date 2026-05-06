from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        substituent_atoms = {
            'F': 9,
            'Cl': 17,
            'O': 8,
            'N': 7,
            'C': 6,
        }

        original_smi = Chem.MolToSmiles(mol_obj)
        derivatives = set()

        for atom_idx in range(mol_obj.GetNumAtoms()):
            atom = mol_obj.GetAtomWithIdx(atom_idx)
            num_implicit_h = atom.GetNumImplicitHs()
            if num_implicit_h > 0:
                for sub_name, atomic_num in substituent_atoms.items():
                    try:
                        rw_mol = Chem.RWMol(mol_obj)
                        new_atom_idx = rw_mol.AddAtom(Chem.Atom(atomic_num))
                        rw_mol.AddBond(atom_idx, new_atom_idx, Chem.BondType.SINGLE)
                        Chem.SanitizeMol(rw_mol)
                        new_smi = Chem.MolToSmiles(rw_mol)
                        if new_smi and new_smi != original_smi:
                            derivatives.add(new_smi)
                    except Exception:
                        pass

        return sorted(list(derivatives))
    except Exception as e:
        print(e)
        return None
