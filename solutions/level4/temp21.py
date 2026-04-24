from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):
    """给定分子 → 判断是否含酮基 → 若有 → 醇缩合 → 计算分子量。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含酮基
        pattern = Chem.MolFromSmarts('[#6][CX3](=O)[#6]')
        has_ketone = mol_obj.HasSubstructMatch(pattern)

        if not has_ketone:
            return None

        # Step 2: 醇缩合（alpha-碳与酮基缩合形成 beta-羟基酮）
        rxn = AllChem.ReactionFromSmarts('[CH2:1][C:2](=O)[#6:3]>>[CH:1]=[C:2](O)[#6:3]')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        # Step 3: 计算分子量
        mol_wt = rdMolDescriptors.CalcExactMolWt(product)

        return {
            "has_ketone": has_ketone,
            "product": product_smiles,
            "molecular_weight": round(mol_wt, 4)
        }
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CC(=O)CC"
    print(f"result: {level_function(smiles)}")
