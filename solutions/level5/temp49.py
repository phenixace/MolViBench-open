from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors, Lipinski
import numpy as np


def level_function(mols):
    """从一组候选 → 筛选 Rule of Five → 做聚类 → 每类保留一个代表分子。"""
    try:
        # Step 1: 筛选 Rule of Five
        ro5_mols = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Lipinski.NumHDonors(mol)
            hba = Lipinski.NumHAcceptors(mol)
            if mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
                ro5_mols.append({
                    'smiles': Chem.MolToSmiles(mol),
                    'fp': list(fp),
                    'qed': round(Descriptors.qed(mol), 4)
                })

        if len(ro5_mols) < 2:
            return [d['smiles'] for d in ro5_mols]

        n = len(ro5_mols)
        k = min(5, n)

        # Step 2: K-means 聚类
        X = np.array([d['fp'] for d in ro5_mols], dtype=float)
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

        # Step 3: 每类保留 QED 最高的代表分子
        results = []
        for cluster_id in range(k):
            cluster_mols = [ro5_mols[i] for i in range(n) if labels[i] == cluster_id]
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
                   "c1ccc(O)cc1", "CCCC", "CCN", "c1ccc(F)cc1",
                   "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"]
    result = level_function(smiles_list)
    if result:
        for r in result:
            print(f"  Cluster {r['cluster']}: {r['smiles']} (QED={r['qed']})")
