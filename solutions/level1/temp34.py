from rdkit import Chem

def level_function(mol1, mol2):
    try:
        mol1 = Chem.MolFromSmiles(mol1)
        if mol1 is None:
            return None
        mol2 = Chem.MolFromSmiles(mol2)
        if mol2 is None:
            return None
        return Chem.MolToSmiles(mol1) == Chem.MolToSmiles(mol2)
    except Exception as e:
        print(e)
        return None
