from rdkit import Chem

def level_function(mol):
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return False
        return True
    except Exception as e:
        print(e)
        return False
