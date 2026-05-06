from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol, filename="output.pdb"):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        mol_obj = Chem.AddHs(mol_obj)
        AllChem.EmbedMolecule(mol_obj, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol_obj)
        Chem.MolToPDBFile(mol_obj, filename)
        return filename
    except Exception as e:
        print(e)
        return None
