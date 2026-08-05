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
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=10, randomSeed=42, maxAttempts=5)
        if len(cids) == 0:
            return None
        for cid in cids:
            AllChem.MMFFOptimizeMolecule(mol, confId=cid, maxIters=100)
        cid_list = list(cids)
        n = len(cid_list)
        matrix = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                rmsd = AllChem.GetConformerRMS(mol, cid_list[i], cid_list[j])
                matrix[i][j] = rmsd
                matrix[j][i] = rmsd
        return matrix
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccccc1CCO'
    result = level_function(smiles)
    if result:
        print(f'Output: {len(result)}{len(result)}')
        for row in result:
            print('Output:', [f'{v:.3f}' for v in row])
