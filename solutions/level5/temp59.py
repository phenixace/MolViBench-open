from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Crippen, FilterCatalog, Lipinski
import numpy as np

def level_function(library_smiles, query_smiles, top_k=10, sim_threshold=0.3):
    try:
        query_mol = Chem.MolFromSmiles(query_smiles)
        if query_mol is None:
            return None
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=2048)

        parsed = []
        for smi in library_smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                parsed.append((Chem.MolToSmiles(mol), mol))

        candidates = []
        for smi, mol in parsed:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            sim = DataStructs.TanimotoSimilarity(query_fp, fp)

            candidates.append({
                "smiles": smi,
                "mol": mol,
                "fp": fp,
                "similarity": sim,
                "MW": Descriptors.MolWt(mol),
                "LogP": Crippen.MolLogP(mol),
                "TPSA": Descriptors.TPSA(mol),
                "HBD": Descriptors.NumHDonors(mol),
                "HBA": Descriptors.NumHAcceptors(mol),
                "RotBonds": Descriptors.NumRotatableBonds(mol),
                "QED": Descriptors.qed(mol),
            })

        n_total = len(candidates)

        candidates = [c for c in candidates if c["similarity"] >= sim_threshold]
        n_after_sim = len(candidates)

        filtered = []
        for c in candidates:
            lipinski_ok = (c["MW"] <= 500 and c["LogP"] <= 5 and c["HBD"] <= 5 and c["HBA"] <= 10)
            veber_ok = (c["TPSA"] <= 140 and c["RotBonds"] <= 10)
            if lipinski_ok and veber_ok:
                filtered.append(c)
        n_after_admet = len(filtered)

        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog.FilterCatalog(params)

        clean = []
        for c in filtered:
            entry = catalog.GetFirstMatch(c["mol"])
            if entry is None:
                clean.append(c)
        n_after_pains = len(clean)

        for c in clean:
            tpsa_norm = min(c["TPSA"], 140) / 140.0
            c["score"] = round(0.5 * c["similarity"] + 0.3 * c["QED"] + 0.2 * (1 - tpsa_norm), 4)

        clean.sort(key=lambda x: x["score"], reverse=True)

        top_results = []
        for c in clean[:top_k]:
            top_results.append({
                "smiles": c["smiles"],
                "similarity": round(c["similarity"], 4),
                "score": c["score"],
                "MW": round(c["MW"], 2),
                "LogP": round(c["LogP"], 4),
                "TPSA": round(c["TPSA"], 2),
                "QED": round(c["QED"], 4),
                "HBD": c["HBD"],
                "HBA": c["HBA"],
            })

        return {
            "query": query_smiles,
            "pipeline_summary": {
                "total_input": n_total,
                "after_similarity_filter": n_after_sim,
                "after_ADMET_filter": n_after_admet,
                "after_PAINS_filter": n_after_pains,
                "top_k_returned": len(top_results),
            },
            "top_hits": top_results,
        }
    except Exception as e:
        print(e)
        return None
