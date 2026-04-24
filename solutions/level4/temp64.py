from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs
from rdkit.Chem.MolStandardize import rdMolStandardize
from sklearn.cluster import KMeans
import numpy as np


def level_function(sdf_content, k=5):
    """从 SDF 文件读取分子 → 标准化（去盐、标准化互变异构体）→ 计算 Morgan 指纹 → K-means 聚类（K=5）→ 每类选 QED 最高的代表分子 → 导出结果。"""
    try:
        supplier = Chem.SDMolSupplier()
        supplier.SetData(sdf_content)

        mols = []
        smiles = []
        for mol in supplier:
            if mol is None:
                continue
            # Standardize
            try:
                chooser = rdMolStandardize.LargestFragmentChooser()
                mol = chooser.choose(mol)
                te = rdMolStandardize.TautomerEnumerator()
                mol = te.Canonicalize(mol)
            except Exception:
                pass
            smi = Chem.MolToSmiles(mol)
            mols.append(mol)
            smiles.append(smi)

        if len(mols) < k:
            return None

        # Morgan fingerprints
        fps = []
        for mol in mols:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fps.append(np.array(fp))
        X = np.array(fps)

        # K-means clustering
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = kmeans.fit_predict(X)

        # Select best QED per cluster
        qed_values = [Descriptors.qed(m) for m in mols]
        representatives = []
        for cluster_id in range(k):
            cluster_indices = [i for i, l in enumerate(labels) if l == cluster_id]
            if not cluster_indices:
                continue
            best_idx = max(cluster_indices, key=lambda i: qed_values[i])
            representatives.append({
                "cluster": cluster_id,
                "smiles": smiles[best_idx],
                "QED": round(qed_values[best_idx], 4),
                "cluster_size": len(cluster_indices)
            })

        return {
            "total_molecules": len(mols),
            "k": k,
            "representatives": representatives
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    # Create test SDF
    test_smiles = ["c1ccccc1", "CCO", "CC(=O)O", "c1ccncc1", "CCCCC",
                   "c1ccc(O)cc1", "c1ccc(N)cc1", "CC(C)C", "CCC(=O)O", "c1ccoc1"]
    sdf_str = ""
    for smi in test_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            AllChem.Compute2DCoords(mol)
            sdf_str += Chem.MolToMolBlock(mol) + "$$$$\n"
    result = level_function(sdf_str, k=3)
    if result:
        print(f"Representatives:")
        for r in result['representatives']:
            print(f"  Cluster {r['cluster']}: {r['smiles']} (QED={r['QED']})")
