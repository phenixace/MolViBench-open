from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


def level_function(mol):
    """给定分子 → 判断是否含芳环 → 若有 → 双硝化 → 计算分子式。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含芳环
        num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol_obj)
        has_aromatic_ring = num_aromatic_rings > 0

        if not has_aromatic_ring:
            return None

        # Step 2: 双硝化 (在芳环上替换两个氢为硝基)
        rxn = AllChem.ReactionFromSmarts('[cH:1]>>[c:1][N+](=O)[O-]')
        products = rxn.RunReactants((mol_obj,))
        if not products or len(products) < 1:
            return None

        # 第一次硝化
        first_product = products[0][0]
        Chem.SanitizeMol(first_product)

        # 第二次硝化
        products2 = rxn.RunReactants((first_product,))
        if not products2:
            return None

        final_product = products2[0][0]
        Chem.SanitizeMol(final_product)
        product_smiles = Chem.MolToSmiles(final_product)

        # Step 3: 计算分子式
        formula = rdMolDescriptors.CalcMolFormula(final_product)

        return {
            "has_aromatic_ring": has_aromatic_ring,
            "product": product_smiles,
            "formula": formula
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccccc1"
    print(f"result: {level_function(smiles)}")
