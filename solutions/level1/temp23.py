from rdkit import Chem

def level_function(mol):



    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return any([atom.GetFormalCharge() != 0 for atom in mol.GetAtoms()])
    except Exception as e:
        print(e)
        return None

if __name__ == '__main__':
    smiles = 'CCOc1ccccc1[O-]'
    print(f'Output: {level_function(smiles)}')
