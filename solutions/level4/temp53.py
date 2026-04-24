from rdkit import Chem
from rdkit.Chem import Descriptors, FilterCatalog


def level_function(smiles_list):
    """给定分子列表 → 逐个检测是否含 PAINS 子结构 → 过滤掉含 PAINS 的分子 → 对剩余分子计算 QED → 返回 QED 最高的前 3 个。"""
    try:
        # Set up PAINS filter
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog.FilterCatalog(params)

        passed = []
        pains_count = 0

        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            entry = catalog.GetFirstMatch(mol)
            if entry is not None:
                pains_count += 1
                continue
            qed = Descriptors.qed(mol)
            passed.append({
                "smiles": Chem.MolToSmiles(mol),
                "QED": round(qed, 4)
            })

        # Sort by QED descending
        passed.sort(key=lambda x: x["QED"], reverse=True)

        return {
            "total_input": len(smiles_list),
            "pains_filtered": pains_count,
            "remaining": len(passed),
            "top3": passed[:3]
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    mols = ["c1ccccc1", "CC(=O)Nc1ccccc1", "CCO", "c1ccc(O)cc1",
            "CC(C)Cc1ccc(C(C)C(=O)O)cc1"]
    result = level_function(mols)
    if result:
        print(f"Filtered: {result['pains_filtered']}, Top-3:")
        for m in result['top3']:
            print(f"  {m['smiles']}: QED={m['QED']}")
