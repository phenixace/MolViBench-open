from rdkit import Chem
from rdkit.Chem import Recap

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        tree = Recap.RecapDecompose(mol_obj)
        leaves = tree.GetLeaves()
        fragments = sorted(list(leaves.keys()))
        return fragments
    except Exception as e:
        print(e)
        return None
