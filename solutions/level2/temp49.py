from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def level_function(mols, k=3, max_iter=100):
    try:
        fps = []
        valid_smiles = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            fps.append(list(fp))
            valid_smiles.append(smi)
        if len(fps) < k:
            return None

        X = np.array(fps, dtype=float)
        n = X.shape[0]

        rng = np.random.RandomState(42)
        indices = rng.choice(n, k, replace=False)
        centroids = X[indices].copy()

        labels = np.zeros(n, dtype=int)
        for _ in range(max_iter):
            dists = np.linalg.norm(X[:, None, :] - centroids[None, :, :], axis=2)
            new_labels = np.argmin(dists, axis=1)
            if np.array_equal(new_labels, labels):
                break
            labels = new_labels
            for j in range(k):
                members = X[labels == j]
                if len(members) > 0:
                    centroids[j] = members.mean(axis=0)

        clusters = {}
        for i, label in enumerate(labels):
            clusters.setdefault(int(label), []).append(valid_smiles[i])
        return clusters
    except Exception as e:
        print(e)
        return None
