from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol):
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        if result != 0:
            return None
        AllChem.MMFFOptimizeMolecule(mol)
        conf = mol.GetConformer()
        coords = []
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            pos = conf.GetAtomPosition(i)
            coords.append((atom.GetSymbol(), pos.x, pos.y, pos.z))
        return coords
    except Exception as e:
        print(e)
        return None
