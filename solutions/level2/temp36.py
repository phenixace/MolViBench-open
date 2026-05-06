from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol, threshold=1000.0):
    try:
        mol = Chem.MolFromSmiles(mol)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        mp = AllChem.MMFFGetMoleculeProperties(mol)
        ff = AllChem.MMFFGetMoleculeForceField(mol, mp)
        if ff is None:
            ff = AllChem.UFFGetMoleculeForceField(mol)
        energy = ff.CalcEnergy()
        return energy < threshold
    except Exception as e:
        print(e)
        return None
