from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, FilterCatalog

BRENK_SMARTS = {
    "aldehyde": "[CH1](=O)", "epoxide": "C1OC1", "peroxide": "OO",
    "azide": "N=[N+]=[N-]", "disulfide": "SS", "nitro": "[N+](=O)[O-]",
}

def level_function(query_smiles, library_smiles, top_k=10):
    try:
        query_mol = Chem.MolFromSmiles(query_smiles)
        if query_mol is None:
            return None

        query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=2048)

        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog.FilterCatalog(params)

        hits = []
        for smi in library_smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            sim = DataStructs.TanimotoSimilarity(query_fp, fp)
            hits.append((Chem.MolToSmiles(mol), mol, sim))

        hits.sort(key=lambda x: x[2], reverse=True)
        top_hits = hits[:top_k]

        clean_results = []
        for smi, mol, sim in top_hits:
            if catalog.GetFirstMatch(mol) is not None:
                continue

            brenk_pass = True
            for name, smarts in BRENK_SMARTS.items():
                pattern = Chem.MolFromSmarts(smarts)
                if pattern and mol.HasSubstructMatch(pattern):
                    brenk_pass = False
                    break
            if not brenk_pass:
                continue

            clean_results.append({
                "smiles": smi,
                "similarity": round(sim, 4)
            })

        return {
            "query": Chem.MolToSmiles(query_mol),
            "top_k_hits": len(top_hits),
            "after_filtering": len(clean_results),
            "clean_molecules": clean_results
        }
    except Exception as e:
        print(e)
        return None
