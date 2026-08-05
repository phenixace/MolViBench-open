from rdkit import Chem

def level_function(mol):



    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        pattern = Chem.MolFromSmarts("n")
        return mol.HasSubstructMatch(pattern)
    except Exception as e:
        print(e)
        return None

if __name__ == '__main__':
    smiles = 'C(N)C[C@H](F)C(=O)O'
    print(f'Output: {level_function(smiles)}')
