from rdkit import Chem
from rdkit.Chem import FilterCatalog

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
        catalog = FilterCatalog.FilterCatalog(params)

        entries = catalog.GetMatches(mol_obj)
        alerts = []
        for entry in entries:
            alerts.append({
                "description": entry.GetDescription(),
            })

        mlsmr_patterns = {
            "alkyl_halide": "[CX4][Cl,Br,I]",
            "acyl_halide": "C(=O)[Cl,Br,I]",
            "sulfonyl_halide": "S(=O)(=O)[Cl,Br,I]",
            "anhydride": "C(=O)OC(=O)",
            "aziridine": "C1NC1",
            "thiol": "[SH]",
            "michael_acceptor_ester": "C=CC(=O)O",
        }

        for name, smarts in mlsmr_patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol_obj.HasSubstructMatch(pattern):
                alerts.append({"description": f"MLSMR_{name}"})

        return {
            "num_alerts": len(alerts),
            "alerts": [a["description"] for a in alerts],
            "passes_MLSMR": len(alerts) == 0
        }
    except Exception as e:
        print(e)
        return None
