from rdkit import Chem

def level_function(mol):



    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        atoms = mol.GetAtoms()
        return sum(1 for atom in atoms if atom.GetSymbol() == "S")
    except Exception as e:
        print(e)
        return None

if __name__ == '__main__':
    smiles = 'SC(N)C[C@H](F)C(=O)O'
    print(f'Output: {level_function(smiles)}')
