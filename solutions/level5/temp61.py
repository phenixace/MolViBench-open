import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Crippen, FilterCatalog
from rdkit.ML.Cluster import Butina

def level_function(query_smiles, library_smiles, top_n=100, cluster_cutoff=0.4):
    try:
        query_mol = Chem.MolFromSmiles(query_smiles)
        if query_mol is None:
            return None
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=2048)

        candidates = []
        for smi in library_smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            sim = DataStructs.TanimotoSimilarity(query_fp, fp)
            candidates.append({
                "smiles": Chem.MolToSmiles(mol),
                "mol": mol,
                "fp": fp,
                "similarity": sim,
            })

        n_parsed = len(candidates)

        candidates.sort(key=lambda x: x["similarity"], reverse=True)
        candidates = candidates[:top_n]
        n_top = len(candidates)

        admet_pass = []
        for c in candidates:
            mol = c["mol"]
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            tpsa = Descriptors.TPSA(mol)
            rotb = Descriptors.NumRotatableBonds(mol)

            lipinski = (mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10)
            veber = (tpsa <= 140 and rotb <= 10)

            if lipinski and veber:
                c["MW"] = mw
                c["LogP"] = logp
                c["TPSA"] = tpsa
                c["QED"] = Descriptors.qed(mol)
                admet_pass.append(c)
        n_admet = len(admet_pass)

        pains_params = FilterCatalog.FilterCatalogParams()
        pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        pains_catalog = FilterCatalog.FilterCatalog(pains_params)

        clean = []
        for c in admet_pass:
            if pains_catalog.GetFirstMatch(c["mol"]) is None:
                clean.append(c)
        n_clean = len(clean)

        if not clean:
            return {
                "query": query_smiles,
                "n_library": n_parsed,
                "n_top_similar": n_top,
                "n_after_ADMET": n_admet,
                "n_after_PAINS": n_clean,
                "n_clusters": 0,
                "representatives": [],
            }

        fps = [c["fp"] for c in clean]
        n = len(fps)
        dists = []
        for i in range(1, n):
            for j in range(i):
                dists.append(1.0 - DataStructs.TanimotoSimilarity(fps[i], fps[j]))

        clusters = Butina.ClusterData(dists, n, cluster_cutoff, isDistData=True)

        representatives = []
        for cluster in clusters:
            cluster_mols = [clean[i] for i in cluster]
            cluster_mols.sort(key=lambda x: x["similarity"], reverse=True)
            best = cluster_mols[0]
            representatives.append({
                "smiles": best["smiles"],
                "similarity": round(best["similarity"], 4),
                "MW": round(best["MW"], 2),
                "LogP": round(best["LogP"], 4),
                "TPSA": round(best["TPSA"], 2),
                "QED": round(best["QED"], 4),
                "cluster_size": len(cluster),
            })

        representatives.sort(key=lambda x: x["similarity"], reverse=True)

        return {
            "query": query_smiles,
            "pipeline_summary": {
                "n_library": n_parsed,
                "n_top_similar": n_top,
                "n_after_ADMET": n_admet,
                "n_after_PAINS": n_clean,
                "n_clusters": len(clusters),
                "n_representatives": len(representatives),
            },
            "representatives": representatives,
        }
    except Exception as e:
        print(e)
        return None
