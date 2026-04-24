from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem import RWMol


def level_function(mol):
    """给定一个药物候选分子，生成环开/环合异构体。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        original_smi = Chem.MolToSmiles(mol_obj)
        isomers = set()

        ring_info = mol_obj.GetRingInfo()
        bond_rings = ring_info.BondRings()

        # Strategy 1: Ring opening - break one bond in each ring
        for ring_bonds in bond_rings:
            for bond_idx in ring_bonds:
                try:
                    rw_mol = Chem.RWMol(mol_obj)
                    bond = rw_mol.GetBondWithIdx(bond_idx)
                    begin_idx = bond.GetBeginAtomIdx()
                    end_idx = bond.GetEndAtomIdx()

                    # Only break single bonds
                    if bond.GetBondType() == Chem.BondType.SINGLE:
                        rw_mol.RemoveBond(begin_idx, end_idx)
                        # Add hydrogens to satisfy valence
                        try:
                            Chem.SanitizeMol(rw_mol)
                            new_smi = Chem.MolToSmiles(rw_mol)
                            if new_smi and new_smi != original_smi:
                                # Verify it's a valid molecule
                                check = Chem.MolFromSmiles(new_smi)
                                if check is not None:
                                    isomers.add(new_smi)
                        except Exception:
                            pass
                except Exception:
                    pass

        # Strategy 2: Ring closure - form new bonds between atoms
        # that are 3-5 bonds apart (to form 4-6 membered rings)
        try:
            dist_matrix = Chem.GetDistanceMatrix(mol_obj)
            num_atoms = mol_obj.GetNumAtoms()

            for i in range(num_atoms):
                for j in range(i + 1, num_atoms):
                    dist = int(dist_matrix[i][j])
                    if 3 <= dist <= 5:
                        atom_i = mol_obj.GetAtomWithIdx(i)
                        atom_j = mol_obj.GetAtomWithIdx(j)

                        # Check if atoms have implicit hydrogens
                        if atom_i.GetNumImplicitHs() > 0 and atom_j.GetNumImplicitHs() > 0:
                            # Check no existing bond
                            if mol_obj.GetBondBetweenAtoms(i, j) is None:
                                try:
                                    rw_mol = Chem.RWMol(mol_obj)
                                    rw_mol.AddBond(i, j, Chem.BondType.SINGLE)
                                    Chem.SanitizeMol(rw_mol)
                                    new_smi = Chem.MolToSmiles(rw_mol)
                                    if new_smi and new_smi != original_smi:
                                        check = Chem.MolFromSmiles(new_smi)
                                        if check is not None:
                                            isomers.add(new_smi)
                                except Exception:
                                    pass
        except Exception:
            pass

        # Strategy 3: Replace ring atoms with different heteroatoms
        atom_rings = ring_info.AtomRings()
        for ring_atoms in atom_rings:
            for atom_idx in ring_atoms:
                atom = mol_obj.GetAtomWithIdx(atom_idx)
                orig_num = atom.GetAtomicNum()
                for new_num in [6, 7, 8, 16]:
                    if new_num != orig_num:
                        try:
                            rw_mol = Chem.RWMol(mol_obj)
                            rw_mol.GetAtomWithIdx(atom_idx).SetAtomicNum(new_num)
                            Chem.SanitizeMol(rw_mol)
                            new_smi = Chem.MolToSmiles(rw_mol)
                            if new_smi and new_smi != original_smi:
                                check = Chem.MolFromSmiles(new_smi)
                                if check is not None:
                                    isomers.add(new_smi)
                        except Exception:
                            pass

        return sorted(list(isomers))
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    result = level_function("C1CCCCC1CC")
    print(f"result: {result}")
