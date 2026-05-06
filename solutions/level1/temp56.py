from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        scaffold = MurckoScaffold.GetScaffoldForMol(mol_obj)
        return Chem.MolToSmiles(scaffold)
    except Exception as e:
        print(e)
        return None
