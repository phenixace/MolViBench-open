from rdkit import Chem

def level_function(mol1, mol2):

    try:
        m1 = Chem.MolFromSmiles(mol1)
        m2 = Chem.MolFromSmiles(mol2)
        if m1 is None or m2 is None: return False



        s1_no = Chem.MolToSmiles(m1, isomericSmiles=False)
        s2_no = Chem.MolToSmiles(m2, isomericSmiles=False)
        if s1_no != s2_no:
            return False


        s1_stereo = Chem.MolToSmiles(m1, isomericSmiles=True)
        s2_stereo = Chem.MolToSmiles(m2, isomericSmiles=True)
        if s1_stereo == s2_stereo:
            return False



        m1_inverted = Chem.Mol(m1)
        centers = Chem.FindMolChiralCenters(m1_inverted, includeUnassigned=True)


        if len(centers) < 1:
            return False

        for idx, _ in centers:
            atom = m1_inverted.GetAtomWithIdx(idx)
            tag = atom.GetChiralTag()
            if tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
            elif tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)


        s1_mirror = Chem.MolToSmiles(m1_inverted, isomericSmiles=True)

        if s1_mirror == s2_stereo:
            return False



        return True

    except Exception as e:
        print(f"Error: {e}")
        return False

if __name__ == '__main__':
    mol1 = '[C@@H](F)(C)[C@@H](O)C'
    mol2 = '[C@@H](F)(C)[C@H](O)C'
    print(f'Output: {level_function(mol1, mol2)}')
