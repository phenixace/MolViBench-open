from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(filename):
    try:
        mol = Chem.MolFromPDBFile(filename, removeHs=True)
        if mol is None:
            return None
        smiles = Chem.MolToSmiles(mol)
        return smiles
    except Exception as e:
        print(e)
        return None
