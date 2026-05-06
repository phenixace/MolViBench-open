from rdkit import Chem
from rdkit.Chem import Descriptors

def level_function(mol):
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return Descriptors.qed(mol)
    except Exception as e:
        print(e)
        return None
