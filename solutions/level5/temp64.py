import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem


def level_function(train_smiles, new_smiles, fp_radius=2, fp_bits=2048, threshold_percentile=95):

    try:

        train_mols = [Chem.MolFromSmiles(s) for s in train_smiles]
        train_mols = [m for m in train_mols if m is not None]
        if len(train_mols) < 3:
            return None

        train_fps = [AllChem.GetMorganFingerprintAsBitVect(m, fp_radius, nBits=fp_bits) for m in train_mols]



        train_nn_dists = []
        for i in range(len(train_fps)):
            min_dist = float('inf')
            for j in range(len(train_fps)):
                if i == j:
                    continue
                d = 1.0 - DataStructs.TanimotoSimilarity(train_fps[i], train_fps[j])
                if d < min_dist:
                    min_dist = d
            train_nn_dists.append(min_dist)

        train_nn_dists = np.array(train_nn_dists)



        mean_dist = float(np.mean(train_nn_dists))
        std_dist = float(np.std(train_nn_dists))
        threshold = float(np.percentile(train_nn_dists, threshold_percentile))


        train_fp_arrays = np.array([list(fp) for fp in train_fps], dtype=float)
        centroid = np.mean(train_fp_arrays, axis=0)

        def centroid_distance(fp_array):
            return float(np.sqrt(np.sum((fp_array - centroid) ** 2)))

        train_centroid_dists = [centroid_distance(fp_arr) for fp_arr in train_fp_arrays]
        centroid_threshold = float(np.percentile(train_centroid_dists, threshold_percentile))


        results = []
        for smi in new_smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                results.append({
                    "smiles": smi,
                    "valid": False,
                    "in_AD": False,
                    "reason": "Invalid SMILES",
                })
                continue

            fp = AllChem.GetMorganFingerprintAsBitVect(mol, fp_radius, nBits=fp_bits)


            nn_dist = min(1.0 - DataStructs.TanimotoSimilarity(fp, tfp) for tfp in train_fps)


            fp_arr = np.array(list(fp), dtype=float)
            c_dist = centroid_distance(fp_arr)


            nn_in_ad = nn_dist <= threshold
            centroid_in_ad = c_dist <= centroid_threshold
            in_ad = nn_in_ad and centroid_in_ad

            results.append({
                "smiles": smi,
                "valid": True,
                "nn_distance": round(nn_dist, 4),
                "nn_threshold": round(threshold, 4),
                "nn_in_AD": nn_in_ad,
                "centroid_distance": round(c_dist, 4),
                "centroid_threshold": round(centroid_threshold, 4),
                "centroid_in_AD": centroid_in_ad,
                "in_AD": in_ad,
            })

        n_in_ad = sum(1 for r in results if r.get("in_AD", False))

        return {
            "training_set": {
                "n_molecules": len(train_mols),
                "nn_distance_mean": round(mean_dist, 4),
                "nn_distance_std": round(std_dist, 4),
                "nn_threshold": round(threshold, 4),
                "centroid_threshold": round(centroid_threshold, 4),
            },
            "evaluation": {
                "n_query": len(new_smiles),
                "n_in_AD": n_in_ad,
                "n_outside_AD": len(results) - n_in_ad,
                "fraction_in_AD": round(n_in_ad / len(results), 4) if results else 0.0,
            },
            "molecule_results": results,
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    train = ['c1ccccc1', 'c1ccc(O)cc1', 'c1ccc(N)cc1', 'c1ccc(F)cc1', 'c1ccc(Cl)cc1', 'c1ccc(Br)cc1', 'c1ccc(C)cc1', 'c1ccc(OC)cc1', 'c1ccc(C(=O)O)cc1', 'c1ccc(C#N)cc1']
    new = ['c1ccc(CC)cc1', 'c1ccncc1', 'C1CC2CCC3CCCC4CCCCC4C3C2C1', 'CCCCCCCCCCCCCC']
    result = level_function(train, new)
    if result:
        print(f"Output: {result['training_set']}")
        print(f"Output: {result['evaluation']}")
        for r in result['molecule_results']:
            print(f"Output: {r['smiles']}{r.get('in_AD')}{r.get('nn_distance')}")
