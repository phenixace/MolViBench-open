from rdkit import Chem

def level_function(mol):



    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return mol.GetRingInfo().NumRings()
    except Exception as e:
        print(e)
        return None

if __name__ == '__main__':
    smiles = 'CCO'
    print(f'Output: {level_function(smiles)}')
