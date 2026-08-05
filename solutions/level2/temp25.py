from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol1, mol2):

    try:
        m1 = Chem.MolFromSmiles(mol1)
        m2 = Chem.MolFromSmiles(mol2)
        if m1 is None or m2 is None:
            return False

        if Chem.MolToSmiles(m1) != Chem.MolToSmiles(m2):
            return False
        m1 = Chem.AddHs(m1)
        m2 = Chem.AddHs(m2)
        if m1.GetNumHeavyAtoms() > 50 or m2.GetNumHeavyAtoms() > 50:
            return False
        res1 = AllChem.EmbedMolecule(m1, randomSeed=42, maxAttempts=5)
        if res1 == -1:
            return False
        res2 = AllChem.EmbedMolecule(m2, randomSeed=42, maxAttempts=5)
        if res2 == -1:
            return False
        AllChem.MMFFOptimizeMolecule(m1, maxIters=100)
        AllChem.MMFFOptimizeMolecule(m2, maxIters=100)

        combined = Chem.RWMol(m1)
        conf2 = m2.GetConformer()
        conf1 = combined.GetConformer()
        rmsd = AllChem.GetBestRMS(m1, m2)
        return rmsd < 0.5
    except Exception as e:
        print(e)
        return False


if __name__ == '__main__':
    smiles1 = 'CCO'
    smiles2 = 'OCCC'
    result = level_function(smiles1, smiles2)
    print(f'Output: {result}')
