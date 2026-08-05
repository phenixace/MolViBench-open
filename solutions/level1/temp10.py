from rdkit import Chem

def level_function(mol):

    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return mol.HasSubstructMatch(Chem.MolFromSmarts("[F,Cl,Br,I]"))
    except Exception as e:
        print(e)
        return None

if __name__ == '__main__':
    smiles = 'CC(F)C(Cl)Br'
    print(f'Output: {level_function(smiles)}')
