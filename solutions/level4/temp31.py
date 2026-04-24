from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


def level_function(mol):
    """给定分子 → 判断是否含苯酚 → 若有 → 磺化 → 计算分子式。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含苯酚
        pattern = Chem.MolFromSmarts("[OH]c1ccccc1")
        has_phenol = mol_obj.HasSubstructMatch(pattern)

        if not has_phenol:
            return None

        # Step 2: 磺化 (芳环上添加 -SO3H)
        rxn = AllChem.ReactionFromSmarts('[cH:1]>>[c:1]S(=O)(=O)O')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        # Step 3: 计算分子式
        formula = rdMolDescriptors.CalcMolFormula(product)

        return {
            "has_phenol": has_phenol,
            "product": product_smiles,
            "formula": formula
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "Oc1ccccc1"
    print(f"result: {level_function(smiles)}")
