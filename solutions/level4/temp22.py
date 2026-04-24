from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):
    """给定分子 → 判断是否含羟基 → 若有 → 醚化 → 计算 TPSA。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含羟基
        pattern = Chem.MolFromSmarts('[OX2H]')
        has_hydroxyl = mol_obj.HasSubstructMatch(pattern)

        if not has_hydroxyl:
            return None

        # Step 2: 醚化（将 -OH 转化为 -OCH3）
        rxn = AllChem.ReactionFromSmarts('[C:1][OH]>>[C:1]OC')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        # Step 3: 计算 TPSA
        tpsa = rdMolDescriptors.CalcTPSA(product)

        return {
            "has_hydroxyl": has_hydroxyl,
            "product": product_smiles,
            "tpsa": round(tpsa, 4)
        }
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CCO"
    print(f"result: {level_function(smiles)}")
