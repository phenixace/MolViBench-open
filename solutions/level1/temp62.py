from rdkit import Chem

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        return Chem.MolToSmarts(mol_obj)
    except Exception as e:
        print(e)
        return None
