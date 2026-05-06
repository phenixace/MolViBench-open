from rdkit import Chem

def level_function(mol, filename="output.mol"):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        Chem.MolToMolFile(mol_obj, filename)
        return filename
    except Exception as e:
        print(e)
        return None
