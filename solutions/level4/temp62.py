from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdFMCS


def level_function(smiles_list):
    """给定一组分子 → 两两计算相似度 → 构建相似度图 → 找出相似度最高的分子对 → 计算该对的 MCS。"""
    try:
        mols = []
        fps = []
        valid_smiles = []
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                mols.append(mol)
                fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048))
                valid_smiles.append(Chem.MolToSmiles(mol))

        if len(mols) < 2:
            return None

        # Pairwise similarity
        n = len(mols)
        best_sim = -1
        best_pair = (0, 1)
        sim_matrix = []
        for i in range(n):
            row = []
            for j in range(n):
                if i == j:
                    row.append(1.0)
                elif j > i:
                    sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
                    row.append(round(sim, 4))
                    if sim > best_sim:
                        best_sim = sim
                        best_pair = (i, j)
                else:
                    row.append(sim_matrix[j][i])
            sim_matrix.append(row)

        # Compute MCS for most similar pair
        i, j = best_pair
        mcs_result = rdFMCS.FindMCS([mols[i], mols[j]], timeout=10)

        return {
            "num_molecules": n,
            "most_similar_pair": {
                "mol1": valid_smiles[i],
                "mol2": valid_smiles[j],
                "similarity": round(best_sim, 4)
            },
            "MCS": {
                "smarts": mcs_result.smartsString,
                "numAtoms": mcs_result.numAtoms,
                "numBonds": mcs_result.numBonds
            },
            "similarity_matrix": sim_matrix
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    mols = ["c1ccccc1O", "c1ccccc1N", "CCCCCC", "c1ccncc1", "c1ccc(O)c(O)c1"]
    result = level_function(mols)
    if result:
        p = result['most_similar_pair']
        print(f"Most similar: {p['mol1']} & {p['mol2']}, sim={p['similarity']}")
        print(f"MCS: {result['MCS']['smarts']}")
