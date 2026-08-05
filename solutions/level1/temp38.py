from rdkit import Chem

def level_function(mol, substructure):



    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        substructure = Chem.MolFromSmarts(substructure)

        if substructure is None:
            return None

        return mol.HasSubstructMatch(substructure)
    except Exception as e:
        print(e)
        return None

if __name__ == '__main__':
    smiles = 'C(N)C[C@H](F)C(=O)O'
    print(f"Output: {level_function(smiles, '[NX3;H2]')}")
