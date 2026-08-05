from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np


def level_function(mols):

    try:
        fps = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            fps.append(list(fp))
        if len(fps) < 2:
            return None
        X = np.array(fps, dtype=float)

        X_centered = X - X.mean(axis=0)
        cov = np.cov(X_centered, rowvar=False)
        eigenvalues, eigenvectors = np.linalg.eigh(cov)

        idx = np.argsort(eigenvalues)[::-1]
        eigenvectors = eigenvectors[:, idx]

        X_pca = X_centered @ eigenvectors[:, :2]
        return X_pca.tolist()
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles_list = ['CCO', 'c1ccccc1', 'CC(=O)O', 'CCCC', 'c1ccc(O)cc1']
    result = level_function(smiles_list)
    print(f'Output: {result}')
