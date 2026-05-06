from rdkit import Chem

def level_function(mol):
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return mol.HasSubstructMatch(Chem.MolFromSmarts("[NX3;H2]"))
    except Exception as e:
        print(e)
        return None
    
