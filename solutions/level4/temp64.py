from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs
from rdkit.Chem.MolStandardize import rdMolStandardize
from sklearn.cluster import KMeans
import numpy as np

def level_function(sdf_content, k=5):
    try:
        supplier = Chem.SDMolSupplier()
        supplier.SetData(sdf_content)

        mols = []
        smiles = []
        for mol in supplier:
            if mol is None:
                continue
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

        fps = []
        for mol in mols:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fps.append(np.array(fp))
        X = np.array(fps)

        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = kmeans.fit_predict(X)

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
