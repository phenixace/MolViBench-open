from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors
import numpy as np


def level_function(mols, k=3):
    """给定一组分子，进行聚类并选择每类中 QED 最高的分子。"""
    try:
        mol_data = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            qed = Descriptors.qed(mol)
            mol_data.append({
                'smiles': Chem.MolToSmiles(mol),
                'fp': list(fp),
                'qed': round(qed, 4)
            })

        n = len(mol_data)
        if n < k:
            k = n

        X = np.array([d['fp'] for d in mol_data], dtype=float)

        # K-means 聚类
        rng = np.random.RandomState(42)
        indices = rng.choice(n, k, replace=False)
        centroids = X[indices].copy()
        labels = np.zeros(n, dtype=int)

        for _ in range(100):
            dists = np.linalg.norm(X[:, None, :] - centroids[None, :, :], axis=2)
            new_labels = np.argmin(dists, axis=1)
            if np.array_equal(new_labels, labels):
                break
            labels = new_labels
            for j in range(k):
                members = X[labels == j]
                if len(members) > 0:
                    centroids[j] = members.mean(axis=0)

        # 每类选 QED 最高的
        results = []
        for cluster_id in range(k):
            cluster_mols = [mol_data[i] for i in range(n) if labels[i] == cluster_id]
            if cluster_mols:
                best = max(cluster_mols, key=lambda x: x['qed'])
                results.append({
                    'cluster': cluster_id,
                    'smiles': best['smiles'],
                    'qed': best['qed']
                })

        return results
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "c1ccncc1",
                   "c1ccc(O)cc1", "CCCC", "CCN", "c1ccc(F)cc1"]
    result = level_function(smiles_list, k=3)
    if result:
        for r in result:
            print(f"  Cluster {r['cluster']}: {r['smiles']} (QED={r['qed']})")
