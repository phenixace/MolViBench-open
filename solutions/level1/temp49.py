from rdkit import Chem

def level_function(mol):



    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        atoms = mol.GetAtoms()
        atom_count = sum(1 for atom in atoms if atom.GetSymbol() == "H")

        return atom_count

    except Exception as e:
        print(e)
        return None

if __name__ == '__main__':
    smiles = 'SC(N)C[C@H](F)C(=O)O'
    print(f'Output: {level_function(smiles)}')
