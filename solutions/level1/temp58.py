from rdkit import Chem
from rdkit.Chem import GraphDescriptors

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        bertz = GraphDescriptors.BertzCT(mol_obj)
        return round(bertz, 4)
    except Exception as e:
        print(e)
        return None
