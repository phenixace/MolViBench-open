from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem import RWMol


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        original_smi = Chem.MolToSmiles(mol_obj)
        isomers = set()

        ring_info = mol_obj.GetRingInfo()
        bond_rings = ring_info.BondRings()


        for ring_bonds in bond_rings:
            for bond_idx in ring_bonds:
                try:
                    rw_mol = Chem.RWMol(mol_obj)
                    bond = rw_mol.GetBondWithIdx(bond_idx)
                    begin_idx = bond.GetBeginAtomIdx()
                    end_idx = bond.GetEndAtomIdx()


                    if bond.GetBondType() == Chem.BondType.SINGLE:
                        rw_mol.RemoveBond(begin_idx, end_idx)

                        try:
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



        try:
            dist_matrix = Chem.GetDistanceMatrix(mol_obj)
            num_atoms = mol_obj.GetNumAtoms()

            for i in range(num_atoms):
                for j in range(i + 1, num_atoms):
                    dist = int(dist_matrix[i][j])
                    if 3 <= dist <= 5:
                        atom_i = mol_obj.GetAtomWithIdx(i)
                        atom_j = mol_obj.GetAtomWithIdx(j)


                        if atom_i.GetNumImplicitHs() > 0 and atom_j.GetNumImplicitHs() > 0:

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


if __name__ == '__main__':
    result = level_function('C1CCCCC1CC')
    print(f'Output: {result}')
