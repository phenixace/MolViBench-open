from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


def level_function(mols):

    try:
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog(params)

        results = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            entry = catalog.GetFirstMatch(mol)
            if entry is None:
                results.append(Chem.MolToSmiles(mol))
        return results
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles_list = ['CCO', 'c1ccccc1', 'c1ccc2c(c1)c(=O)c1ccccc1o2']
    print(f'Output: {level_function(smiles_list)}')
