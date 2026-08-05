from rdkit import Chem

def level_function(mol):



    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        molblock = Chem.MolToMolBlock(mol)
        return molblock
    except Exception as e:
        print(e)
        return None

if __name__ == '__main__':
    smiles = 'CC[C@H](F)C(=O)O'
    print(f'Output: {level_function(smiles)}')
