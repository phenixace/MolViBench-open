from rdkit import Chem

def level_function(mol1, mol2):
    try:
        m1 = Chem.MolFromSmiles(mol1)
        m2 = Chem.MolFromSmiles(mol2)
        if m1 is None or m2 is None:
            return False

        s1 = Chem.MolToSmiles(m1, isomericSmiles=True)
        s2 = Chem.MolToSmiles(m2, isomericSmiles=True)
        if s1 == s2:
            return False


        s1_no = Chem.MolToSmiles(m1, isomericSmiles=False)
        s2_no = Chem.MolToSmiles(m2, isomericSmiles=False)
        if s1_no != s2_no:
            return False

        m1_inverted = Chem.Mol(m1)
        centers = Chem.FindMolChiralCenters(m1_inverted, includeUnassigned=True)

        if not centers:
            return False

        for idx, _ in centers:
            atom = m1_inverted.GetAtomWithIdx(idx)
            curr_tag = atom.GetChiralTag()
            if curr_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
            elif curr_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)

        s1_inverted = Chem.MolToSmiles(m1_inverted, isomericSmiles=True)

        return s1_inverted == s2

    except Exception as e:
        print(f"Error: {e}")
        return False

if __name__ == "__main__":
    mol1 = "C([C@@H](F)Cl)O"
    mol2 = "C([C@H](F)Cl)O"
    print(f"Output: {level_function(mol1, mol2)}")