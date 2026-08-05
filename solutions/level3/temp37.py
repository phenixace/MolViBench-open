from rdkit import Chem

def level_function(mol1, mol2):

    try:
        m1 = Chem.MolFromSmiles(mol1)
        m2 = Chem.MolFromSmiles(mol2)
        if m1 is None or m2 is None: return False


        if Chem.GetFormalCharge(m1) != Chem.GetFormalCharge(m2):
            return False

        def get_standardized_skeleton(mol):

            new_mol = Chem.AddHs(mol)
            rw_mol = Chem.RWMol(new_mol)

            for atom in rw_mol.GetAtoms():
                atom.SetFormalCharge(0)

            for bond in rw_mol.GetBonds():
                bond.SetBondType(Chem.rdchem.BondType.SINGLE)



            return Chem.MolToSmiles(rw_mol, canonical=True)

        skeleton1 = get_standardized_skeleton(m1)
        skeleton2 = get_standardized_skeleton(m2)

        return skeleton1 == skeleton2

    except Exception as e:
        print(f"Error: {e}")
        return False

if __name__ == '__main__':
    print(f"Output: {level_function('C=C[O-]', '[CH2-]C=O')}")
    print(f"Output: {level_function('CC(C)=O', 'CC(C)=CO')}")
