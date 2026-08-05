from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol, filename="output.pdb"):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        mol_obj = Chem.AddHs(mol_obj)
        if mol_obj.GetNumHeavyAtoms() > 50:
            return None
        res = AllChem.EmbedMolecule(mol_obj, randomSeed=42, maxAttempts=5)
        if res == -1:
            return None
        AllChem.MMFFOptimizeMolecule(mol_obj, maxIters=100)
        Chem.MolToPDBFile(mol_obj, filename)
        return filename
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'CCO'
    result = level_function(smiles, 'output.pdb')
    print(f'Output: {result}')
