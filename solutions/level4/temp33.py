from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


def level_function(mol):
    """给定分子 → 判断是否含醇 → 若有 → 脱水成烯烃 → 计算分子量。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含醇 (非芳香 -OH)
        pattern = Chem.MolFromSmarts('[CX4][OX2H]')
        has_alcohol = mol_obj.HasSubstructMatch(pattern)

        if not has_alcohol:
            return None

        # Step 2: 脱水成烯烃
        rxn = AllChem.ReactionFromSmarts('[CH:1][C:2][OH:3]>>[C:1]=[C:2]')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        # Step 3: 计算分子量
        mol_wt = rdMolDescriptors.CalcExactMolWt(product)

        return {
            "has_alcohol": has_alcohol,
            "product": product_smiles,
            "molecular_weight": round(mol_wt, 4)
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "CCCO"
    print(f"result: {level_function(smiles)}")
