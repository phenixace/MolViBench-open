from rdkit import Chem
from rdkit.Chem import FilterCatalog


def level_function(mol):
    """检测分子是否命中 NIH MLSMR 结构警报（扩展版 PAINS）。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Use RDKit's built-in PAINS filter catalog
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

        # Additional NIH MLSMR alerts (beyond PAINS)
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


if __name__ == "__main__":
    smiles = "O=C1C=CC(=O)C=C1"  # Quinone - known PAINS
    result = level_function(smiles)
    if result:
        print(f"MLSMR passes: {result['passes_MLSMR']}")
        print(f"Alerts: {result['alerts']}")
