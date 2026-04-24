from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import numpy as np


def level_function(smiles_list, n_select=20, seed=42):
    """给定一组分子，使用 MaxMin 算法选择结构多样性最大的子集（N=20）。"""
    try:
        np.random.seed(seed)

        mols = []
        fps = []
        valid = []
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                mols.append(mol)
                fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048))
                valid.append(Chem.MolToSmiles(mol))

        if len(mols) <= n_select:
            return {"selected": valid, "num_selected": len(valid)}

        n = len(fps)

        # MaxMin algorithm
        # Start with the first molecule
        selected_indices = [0]
        remaining = set(range(1, n))

        while len(selected_indices) < n_select and remaining:
            best_idx = None
            best_min_dist = -1

            for i in remaining:
                # Minimum distance to any selected molecule
                min_dist = min(
                    1 - DataStructs.TanimotoSimilarity(fps[i], fps[j])
                    for j in selected_indices
                )
                if min_dist > best_min_dist:
                    best_min_dist = min_dist
                    best_idx = i

            if best_idx is not None:
                selected_indices.append(best_idx)
                remaining.remove(best_idx)

        selected_smiles = [valid[i] for i in selected_indices]

        # Calculate diversity metrics
        selected_fps = [fps[i] for i in selected_indices]
        pairwise_dists = []
        for i in range(len(selected_fps)):
            for j in range(i + 1, len(selected_fps)):
                pairwise_dists.append(
                    1 - DataStructs.TanimotoSimilarity(selected_fps[i], selected_fps[j])
                )

        return {
            "total_input": len(mols),
            "num_selected": len(selected_smiles),
            "selected": selected_smiles,
            "avg_pairwise_distance": round(np.mean(pairwise_dists), 4) if pairwise_dists else 0,
            "min_pairwise_distance": round(min(pairwise_dists), 4) if pairwise_dists else 0
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    library = [f"C{'C'*i}O" for i in range(30)] + ["c1ccccc1", "c1ccncc1", "c1ccoc1"]
    result = level_function(library, n_select=5)
    if result:
        print(f"Selected {result['num_selected']} from {result['total_input']}")
        print(f"Avg distance: {result['avg_pairwise_distance']}")
