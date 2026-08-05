from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol):

    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        if mol.GetNumHeavyAtoms() > 50:
            return None
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=2, randomSeed=42, maxAttempts=5)
        if len(cids) < 2:
            return None
        for cid in cids:
            AllChem.MMFFOptimizeMolecule(mol, confId=cid, maxIters=100)
        rmsd = AllChem.GetConformerRMS(mol, cids[0], cids[1])
        return rmsd
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccccc1CCO'
    result = level_function(smiles)
    print(f'Output: {result}')
