from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


def level_function(mol):
    """给定分子，预测是否符合 PAINS 过滤规则。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # 使用 RDKit 内置的 PAINS 过滤器
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


if __name__ == "__main__":
    smiles1 = "CCO"
    smiles2 = "c1ccc2c(c1)c(=O)c1ccccc1o2"
    print(f"乙醇 PAINS 检查: {level_function(smiles1)}")
    print(f"黄酮 PAINS 检查: {level_function(smiles2)}")
