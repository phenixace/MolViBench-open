from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol):
    """给定一个 SMILES，使用枚举替代法生成所有卤素取代体。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        original_smi = Chem.MolToSmiles(mol_obj)
        halogens = [9, 17, 35, 53]  # F, Cl, Br, I
        halogen_names = {9: 'F', 17: 'Cl', 35: 'Br', 53: 'I'}

        derivatives = set()

        # For each atom with implicit H, replace one H with each halogen
        for atom_idx in range(mol_obj.GetNumAtoms()):
            atom = mol_obj.GetAtomWithIdx(atom_idx)
            if atom.GetNumImplicitHs() > 0:
                for halogen in halogens:
                    try:
                        rw_mol = Chem.RWMol(mol_obj)
                        new_idx = rw_mol.AddAtom(Chem.Atom(halogen))
                        rw_mol.AddBond(atom_idx, new_idx, Chem.BondType.SINGLE)
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


if __name__ == "__main__":
    result = level_function("c1ccccc1")
    print(f"result: {result}")
