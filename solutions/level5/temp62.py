import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Crippen


def level_function(generated_smiles, reference_smiles):

    try:
        n_total = len(generated_smiles)
        if n_total == 0:
            return None


        valid_mols = []
        valid_smiles = []
        for smi in generated_smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                canon = Chem.MolToSmiles(mol)
                valid_mols.append(mol)
                valid_smiles.append(canon)

        validity = len(valid_smiles) / n_total


        unique_smiles = set(valid_smiles)
        uniqueness = len(unique_smiles) / len(valid_smiles) if valid_smiles else 0.0


        ref_canonical = set()
        for smi in reference_smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                ref_canonical.add(Chem.MolToSmiles(mol))

        novel_smiles = unique_smiles - ref_canonical
        novelty = len(novel_smiles) / len(unique_smiles) if unique_smiles else 0.0


        unique_list = list(unique_smiles)
        unique_mols = [Chem.MolFromSmiles(s) for s in unique_list]
        fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in unique_mols if m is not None]

        if len(fps) >= 2:
            dists = []
            for i in range(len(fps)):
                for j in range(i + 1, len(fps)):
                    dists.append(1.0 - DataStructs.TanimotoSimilarity(fps[i], fps[j]))
            internal_diversity = float(np.mean(dists))
        else:
            internal_diversity = 0.0


        mw_vals, logp_vals, tpsa_vals, qed_vals = [], [], [], []
        for mol in unique_mols:
            if mol is not None:
                mw_vals.append(Descriptors.MolWt(mol))
                logp_vals.append(Crippen.MolLogP(mol))
                tpsa_vals.append(Descriptors.TPSA(mol))
                qed_vals.append(Descriptors.qed(mol))

        def stats(vals):
            if not vals:
                return {"mean": 0.0, "std": 0.0, "min": 0.0, "max": 0.0, "median": 0.0}
            arr = np.array(vals)
            return {
                "mean": round(float(np.mean(arr)), 4),
                "std": round(float(np.std(arr)), 4),
                "min": round(float(np.min(arr)), 4),
                "max": round(float(np.max(arr)), 4),
                "median": round(float(np.median(arr)), 4),
            }


        ref_mols = [Chem.MolFromSmiles(s) for s in reference_smiles]
        ref_mols = [m for m in ref_mols if m is not None]
        ref_mw = [Descriptors.MolWt(m) for m in ref_mols]
        ref_logp = [Crippen.MolLogP(m) for m in ref_mols]


        def hist_overlap(vals1, vals2, bins=20):
            if not vals1 or not vals2:
                return 0.0
            all_vals = vals1 + vals2
            lo, hi = min(all_vals), max(all_vals)
            if lo == hi:
                return 1.0
            h1, _ = np.histogram(vals1, bins=bins, range=(lo, hi), density=True)
            h2, _ = np.histogram(vals2, bins=bins, range=(lo, hi), density=True)
            h1 = h1 / (h1.sum() + 1e-10)
            h2 = h2 / (h2.sum() + 1e-10)
            return float(np.sum(np.minimum(h1, h2)))

        mw_overlap = hist_overlap(mw_vals, ref_mw)
        logp_overlap = hist_overlap(logp_vals, ref_logp)

        return {
            "n_generated": n_total,
            "n_valid": len(valid_smiles),
            "n_unique": len(unique_smiles),
            "n_novel": len(novel_smiles),
            "metrics": {
                "validity": round(validity, 4),
                "uniqueness": round(uniqueness, 4),
                "novelty": round(novelty, 4),
                "internal_diversity": round(internal_diversity, 4),
            },
            "property_distributions": {
                "MW": stats(mw_vals),
                "LogP": stats(logp_vals),
                "TPSA": stats(tpsa_vals),
                "QED": stats(qed_vals),
            },
            "distribution_overlap": {
                "MW_overlap": round(mw_overlap, 4),
                "LogP_overlap": round(logp_overlap, 4),
            },
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    generated = ['c1ccccc1', 'c1ccc(O)cc1', 'c1ccc(N)cc1', 'invalid_smiles', 'c1ccc(F)cc1', 'c1ccc(Cl)cc1', 'c1ccccc1', 'CCO', 'CCCO', 'CC(=O)O']
    reference = ['c1ccccc1', 'CCO', 'CC(C)O', 'c1ccncc1']
    result = level_function(generated, reference)
    if result:
        print(f"Output: {result['metrics']['validity']}")
        print(f"Output: {result['metrics']['uniqueness']}")
        print(f"Output: {result['metrics']['novelty']}")
        print(f"Output: {result['metrics']['internal_diversity']}")
