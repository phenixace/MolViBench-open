import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from sklearn.cluster import MiniBatchKMeans
from concurrent.futures import ProcessPoolExecutor
import time

def _process_single_mol(mol_block):
    try:
        mol = Chem.MolFromMolBlock(mol_block)
        if mol is None:
            return None

        chooser = rdMolStandardize.LargestFragmentChooser()
        mol = chooser.choose(mol)

        te = rdMolStandardize.TautomerEnumerator()
        te.SetMaxTautomers(20)
        mol = te.Canonicalize(mol)

        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        arr = np.zeros((2048,), dtype=np.float32)
        DataStructs.ConvertToNumpyArray(fp, arr)

        qed_val = Descriptors.qed(mol)

        return {
            "smiles": Chem.MolToSmiles(mol),
            "fp": arr,
            "qed": qed_val
        }
    except:
        return None

def level_function(sdf_content, k=5, n_workers=4):
    start_time = time.time()

    mol_blocks = [b + "$$$$" for b in sdf_content.split("$$$$\n") if b.strip()]
    total_input = len(mol_blocks)
    print(f"[*] Starting task: Total {total_input} molecule blocks")

    valid_results = []
    print(f"[*] Cleaning structures and extracting features (parallel cores: {n_workers})...")
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        results = list(executor.map(_process_single_mol, mol_blocks))
        valid_results = [r for r in results if r is not None]

    num_valid = len(valid_results)
    if num_valid < k:
        print(f"[Error] Too few valid molecules ({num_valid}), cannot cluster into {k} classes")
        return None

    print(f"[*] Building feature matrix...")
    X = np.ascontiguousarray(
        np.array([item['fp'] for item in valid_results], dtype=np.float32)
    )
    qed_values = np.array([item['qed'] for item in valid_results])
    smiles_list = [item['smiles'] for item in valid_results]

    print(f"[*] Performing clustering calculation (n={num_valid}, dim=2048)...")
    try:
        bsize = max(k * 3, 256)
        kmeans = MiniBatchKMeans(
            n_clusters=k,
            batch_size=min(bsize, num_valid),
            n_init="auto",
            random_state=42,
            max_iter=100
        )
        labels = kmeans.fit_predict(X)
    except Exception as e:
        print(f"[Error] Clustering algorithm crashed: {e}")
        return None

    print(f"[*] Filtering representative molecules for each cluster...")
    representatives = []
    for cluster_id in range(k):
        indices = np.where(labels == cluster_id)[0]
        if len(indices) == 0: continue

        best_idx = indices[np.argmax(qed_values[indices])]

        representatives.append({
            "cluster": int(cluster_id),
            "smiles": smiles_list[best_idx],
            "QED": round(float(qed_values[best_idx]), 4),
            "cluster_size": len(indices)
        })

    duration = time.time() - start_time
    print(f"[*] Task completed! Total time: {duration:.2f} seconds")

    return {
        "stats": {"input": total_input, "valid": num_valid},
        "representatives": representatives
    }

if __name__ == "__main__":
    raw_smiles = ["c1ccccc1", "CCO", "CC(=O)O", "c1ccncc1", "CCCCC",
                  "c1ccc(O)cc1", "c1ccc(N)cc1", "CC(C)C", "CCC(=O)O", "c1ccoc1"]
    test_smiles = raw_smiles

    test_sdf = ""
    for smi in test_smiles:
        m = Chem.MolFromSmiles(smi)
        if m:
            test_sdf += Chem.MolToMolBlock(m) + "$$$$\n"

    res = level_function(test_sdf, k=5, n_workers=4)
    print(f"Output: {res}")
