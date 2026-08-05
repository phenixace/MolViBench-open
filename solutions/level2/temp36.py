from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol, threshold=1000.0):

    try:
        mol = Chem.MolFromSmiles(mol)
        mol = Chem.AddHs(mol)
        if mol.GetNumHeavyAtoms() > 50:
            return None
        res = AllChem.EmbedMolecule(mol, randomSeed=42, maxAttempts=5)
        if res == -1:
            return None
        AllChem.MMFFOptimizeMolecule(mol, maxIters=100)
        mp = AllChem.MMFFGetMoleculeProperties(mol)
        ff = AllChem.MMFFGetMoleculeForceField(mol, mp)
        if ff is None:
            ff = AllChem.UFFGetMoleculeForceField(mol)
        energy = ff.CalcEnergy()
        return energy < threshold
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    result = level_function('CCO', 1000.0)
    print(f'Output: {result}')
