from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol):

    try:
        mol = Chem.MolFromSmiles(mol)
        mol = Chem.AddHs(mol)
        if mol.GetNumHeavyAtoms() > 50:
            return None
        res = AllChem.EmbedMolecule(mol, randomSeed=42, maxAttempts=5)
        if res == -1:
            return None
        AllChem.MMFFOptimizeMolecule(mol, maxIters=100)
        molblock = Chem.MolToMolBlock(mol)
        return molblock
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    result = level_function('c1ccccc1CCO')
    print(f'Output: {result}')
