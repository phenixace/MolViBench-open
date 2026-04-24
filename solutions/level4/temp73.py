from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdFMCS
from rdkit.ML.Cluster import Butina


def level_function(smiles_list, distance_threshold=0.5):
    """给定分子集合 → Butina 聚类 → 对每个聚类提取 MCS → 输出每个 scaffold SMILES 及其在聚类中的出现频率。"""
    try:
        mols = []
        fps = []
        valid = []
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                mols.append(mol)
                fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048))
                valid.append(Chem.MolToSmiles(mol))

        if len(mols) < 2:
            return None

        # Pairwise distances
        n = len(fps)
        dists = []
        for i in range(1, n):
            for j in range(i):
                dists.append(1 - DataStructs.TanimotoSimilarity(fps[i], fps[j]))

        clusters = Butina.ClusterData(dists, n, distance_threshold, isDistData=True)

        result = []
        for cluster_id, cluster_indices in enumerate(clusters):
            cluster_mols = [mols[i] for i in cluster_indices]
            cluster_smiles = [valid[i] for i in cluster_indices]

            # Extract MCS for cluster
            if len(cluster_mols) >= 2:
                mcs = rdFMCS.FindMCS(cluster_mols, timeout=5)
                mcs_smarts = mcs.smartsString
                mcs_atoms = mcs.numAtoms
            else:
                mcs_smarts = Chem.MolToSmarts(cluster_mols[0])
                mcs_atoms = cluster_mols[0].GetNumAtoms()

            result.append({
                "cluster_id": cluster_id,
                "size": len(cluster_indices),
                "mcs_smarts": mcs_smarts,
                "mcs_numAtoms": mcs_atoms,
                "members": cluster_smiles
            })

        return {
            "total_molecules": len(mols),
            "num_clusters": len(clusters),
            "clusters": result
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    mols = ["c1ccccc1O", "c1ccccc1N", "c1ccccc1F",
            "CCCCCC", "CCCCCCC", "c1ccncc1", "c1ccoc1"]
    result = level_function(mols, 0.6)
    if result:
        print(f"Clusters: {result['num_clusters']}")
        for c in result['clusters']:
            print(f"  Cluster {c['cluster_id']}: {c['size']} mols, MCS atoms={c['mcs_numAtoms']}")
