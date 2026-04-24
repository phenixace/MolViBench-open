from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):
    """给定分子 → 判断是否含羧基 → 若有 → 转化为酯 → 计算分子量。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含羧基
        pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H]')
        has_carboxyl = mol_obj.HasSubstructMatch(pattern)

        if not has_carboxyl:
            return None

        # Step 2: 转化为甲酯
        rxn = AllChem.ReactionFromSmarts('[C:1](=O)[OH]>>[C:1](=O)OC')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        # Step 3: 计算分子量
        mol_wt = rdMolDescriptors.CalcExactMolWt(product)

        return {
            "has_carboxyl": has_carboxyl,
            "product": product_smiles,
            "molecular_weight": round(mol_wt, 4)
        }
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CC(=O)O"
    print(f"result: {level_function(smiles)}")
