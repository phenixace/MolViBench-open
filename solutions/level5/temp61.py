import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Crippen, FilterCatalog
from rdkit.ML.Cluster import Butina


def level_function(query_smiles, library_smiles, top_n=100, cluster_cutoff=0.4):
    """相似性搜索驱动的虚拟筛选：给定活性分子 query → Top-100 相似性搜索 → 级联 ADMET 过滤 → Butina 聚类 → 每类选最优 → 输出报告。"""
    try:
        query_mol = Chem.MolFromSmiles(query_smiles)
        if query_mol is None:
            return None
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=2048)

        # Step 1: Parse library and compute similarity
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

        # Step 2: Top-N similarity search
        candidates.sort(key=lambda x: x["similarity"], reverse=True)
        candidates = candidates[:top_n]
        n_top = len(candidates)

        # Step 3: ADMET cascade filtering
        # Lipinski
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

        # PAINS filter
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

        # Step 4: Butina clustering
        fps = [c["fp"] for c in clean]
        n = len(fps)
        dists = []
        for i in range(1, n):
            for j in range(i):
                dists.append(1.0 - DataStructs.TanimotoSimilarity(fps[i], fps[j]))

        clusters = Butina.ClusterData(dists, n, cluster_cutoff, isDistData=True)

        # Step 5: Select best (highest similarity) from each cluster
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


if __name__ == "__main__":
    query = "c1ccc(NC(=O)c2ccccc2)cc1"
    library = [
        "c1ccc(NC(=O)c2ccccc2)cc1", "c1ccc(NC(=O)c2ccncc2)cc1",
        "c1ccc(NC(=O)c2ccc(F)cc2)cc1", "c1ccc(NC(=O)C)cc1",
        "c1ccc(NC(=O)CC)cc1", "c1ccc(O)cc1", "c1ccc(N)cc1",
        "CC(=O)Oc1ccccc1C(=O)O", "c1ccc2c(c1)cc1ccccc12",
        "c1ccc(CC(=O)O)cc1", "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
    ]
    result = level_function(query, library, top_n=10)
    if result:
        print(f"Pipeline: {result['pipeline_summary']}")
        for r in result["representatives"][:5]:
            print(f"  {r['smiles']}: sim={r['similarity']}, QED={r['QED']}")
