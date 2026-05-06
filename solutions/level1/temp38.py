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
