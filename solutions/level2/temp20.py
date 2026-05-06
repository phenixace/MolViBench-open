from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        mol_h = Chem.AddHs(mol_obj)
        result = AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())
        if result == -1:
            return None
        opt_result = AllChem.UFFOptimizeMolecule(mol_h)
        conf = mol_h.GetConformer()
        coords = []
        for i in range(mol_h.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            coords.append((mol_h.GetAtomWithIdx(i).GetSymbol(), pos.x, pos.y, pos.z))
        return coords
    except Exception as e:
        print(e)
        return None
