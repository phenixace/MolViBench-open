from rdkit import Chem
from rdkit.Chem import Crippen

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        contribs = Crippen._GetAtomContribs(mol_obj)
        atom_mr = {i: round(c[1], 4) for i, c in enumerate(contribs)}
        total_mr = round(sum(c[1] for c in contribs), 4)
        return {"atom_contributions": atom_mr, "total_MR": total_mr}
    except Exception as e:
        print(e)
        return None
