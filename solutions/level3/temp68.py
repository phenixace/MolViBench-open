from rdkit import Chem
from rdkit.Chem import FilterCatalog


def get_pains_catalog():
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
    return FilterCatalog.FilterCatalog(params)

PAINS_CATALOG = get_pains_catalog()

def level_function(mol_smiles):
    try:
        mol_obj = Chem.MolFromSmiles(mol_smiles)
        if mol_obj is None: return None


        matches = PAINS_CATALOG.GetMatches(mol_obj)
        alerts = [m.GetDescription() for m in matches]


        mlsmr_patterns = {
            "alkyl_halide": "[CX4][Cl,Br,I]",
            "acyl_halide": "C(=O)[Cl,Br,I]",
            "sulfonyl_halide": "S(=O)(=O)[Cl,Br,I]",
            "anhydride": "C(=O)OC(=O)",
            "aziridine": "C1NC1",
            "thiol": "[SH]",
            "michael_acceptor": "[#6&D2]=[#6&D1]-[C,S,N]=[O,S,N]",
        }

        for name, smarts in mlsmr_patterns.items():
            patt = Chem.MolFromSmarts(smarts)
            if mol_obj.HasSubstructMatch(patt):
                alerts.append(f"MLSMR_{name}")

        return {
            "num_alerts": len(alerts),
            "alerts": list(set(alerts)),
            "passes_MLSMR": len(alerts) == 0
        }
    except Exception as e:
        return None


if __name__ == '__main__':
    smiles = 'O=C1C=CC(=O)C=C1'
    result = level_function(smiles)
    if result:
        print(f"Output: {result['passes_MLSMR']}")
        print(f"Output: {result['alerts']}")
