from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs


def level_function(mols):
    """给定一组分子，基于结构相似性构建一个候选库。"""
    try:
        mol_data = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                canonical = Chem.MolToSmiles(mol)
                mol_data.append({'smiles': canonical, 'fp': fp})

        if len(mol_data) < 2:
            return [d['smiles'] for d in mol_data]

        n = len(mol_data)
        similarity_threshold = 0.4

        # Compute pairwise similarity matrix
        sim_matrix = []
        for i in range(n):
            row = []
            for j in range(n):
                if i == j:
                    row.append(1.0)
                else:
                    sim = DataStructs.TanimotoSimilarity(mol_data[i]['fp'], mol_data[j]['fp'])
                    row.append(sim)
            sim_matrix.append(row)

        # Cluster by similarity: greedy clustering
        clusters = []
        assigned = [False] * n
        for i in range(n):
            if assigned[i]:
                continue
            cluster = [i]
            assigned[i] = True
            for j in range(i + 1, n):
                if not assigned[j]:
                    if sim_matrix[i][j] >= similarity_threshold:
                        cluster.append(j)
                        assigned[j] = True
            clusters.append(cluster)

        # Select the molecule with lowest index as representative from each cluster
        # (could also choose most "central" molecule)
        representatives = []
        for cluster in clusters:
            # Pick the molecule with highest avg similarity to cluster members
            best_idx = cluster[0]
            if len(cluster) > 1:
                best_avg_sim = 0
                for idx in cluster:
                    avg_sim = sum(sim_matrix[idx][j] for j in cluster) / len(cluster)
                    if avg_sim > best_avg_sim:
                        best_avg_sim = avg_sim
                        best_idx = idx
            representatives.append(mol_data[best_idx]['smiles'])

        return representatives
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = [
        "c1ccccc1",
        "c1ccc(C)cc1",
        "c1ccc(O)cc1",
        "CCO",
        "CCCO",
        "CC(=O)O",
        "c1ccncc1",
        "c1ccc(F)cc1",
    ]
    result = level_function(smiles_list)
    print(f"result: {result}")
