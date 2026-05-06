from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def level_function(mol):
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return rdMolDescriptors.CalcNumHBA(mol)
    except Exception as e:
        print(e)
        return None
