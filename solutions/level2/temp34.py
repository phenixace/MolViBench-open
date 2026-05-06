from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol):
    try:
        mol = Chem.MolFromSmiles(mol)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        conf = mol.GetConformer()
        coords = []
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            coords.append((atom.GetSymbol(), round(pos.x, 4), round(pos.y, 4), round(pos.z, 4)))
        return coords
    except Exception as e:
        print(e)
        return None
