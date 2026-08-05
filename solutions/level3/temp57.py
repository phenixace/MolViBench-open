from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS


def level_function(mol1, mol2):

    try:
        m1 = Chem.MolFromSmiles(mol1)
        m2 = Chem.MolFromSmiles(mol2)
        if m1 is None or m2 is None:
            return None


        mcs_result = rdFMCS.FindMCS([m1, m2],
                                     atomCompare=rdFMCS.AtomCompare.CompareElements,
                                     bondCompare=rdFMCS.BondCompare.CompareOrder,
                                     timeout=10)

        if mcs_result.numAtoms == 0:
            return None

        mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
        if mcs_mol is None:
            return None


        match1 = m1.GetSubstructMatch(mcs_mol)
        match2 = m2.GetSubstructMatch(mcs_mol)

        if not match1 or not match2:
            return None


        diff_atoms1 = set(range(m1.GetNumAtoms())) - set(match1)
        diff_atoms2 = set(range(m2.GetNumAtoms())) - set(match2)


        def get_diff_smiles(mol, diff_atoms):
            if not diff_atoms:
                return "[H]"
            try:
                edit = Chem.RWMol(mol)
                atoms_to_remove = sorted(set(range(mol.GetNumAtoms())) - diff_atoms, reverse=True)
                for idx in atoms_to_remove:
                    edit.RemoveAtom(idx)
                try:
                    Chem.SanitizeMol(edit)
                    return Chem.MolToSmiles(edit)
                except Exception:
                    return str(sorted(diff_atoms))
            except Exception:
                return str(sorted(diff_atoms))

        frag1 = get_diff_smiles(m1, diff_atoms1)
        frag2 = get_diff_smiles(m2, diff_atoms2)

        return {
            "core_smarts": mcs_result.smartsString,
            "core_num_atoms": mcs_result.numAtoms,
            "transformation": f"{frag1} >> {frag2}",
            "mol1_diff_atoms": sorted(list(diff_atoms1)),
            "mol2_diff_atoms": sorted(list(diff_atoms2))
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smi1 = 'c1ccc(Cl)cc1'
    smi2 = 'c1ccc(F)cc1'
    result = level_function(smi1, smi2)
    if result:
        print(f"Output: {result['transformation']}")
        print(f"Output: {result['core_smarts']}")
