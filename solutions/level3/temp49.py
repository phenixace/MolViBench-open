from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None


        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog(params)

        entry = catalog.GetFirstMatch(mol_obj)
        if entry is not None:
            return {
                "passes_pains": False,
                "alert_name": entry.GetDescription(),
            }
        else:
            return {
                "passes_pains": True,
                "alert_name": None
            }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles1 = 'CCO'
    smiles2 = 'c1ccc2c(c1)c(=O)c1ccccc1o2'
    print(f'Output: {level_function(smiles1)}')
    print(f'Output: {level_function(smiles2)}')
