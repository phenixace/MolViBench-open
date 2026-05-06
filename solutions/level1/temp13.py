from rdkit import Chem
from rdkit.Chem import rdinchi

def level_function(mol):
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return rdinchi.MolToInchi(mol)[0]
    except Exception as e:
        print(e)
        return None
