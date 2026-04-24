from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):
    """给定分子 → 判断是否含双键 → 若有 → 氢化 → 计算分子量。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含 C=C 双键
        pattern = Chem.MolFromSmarts('[C:1]=[C:2]')
        has_double_bond = mol_obj.HasSubstructMatch(pattern)

        if not has_double_bond:
            return None

        # Step 2: 氢化（将 C=C 还原为 C-C）
        rxn = AllChem.ReactionFromSmarts('[C:1]=[C:2]>>[C:1][C:2]')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        # Step 3: 计算分子量
        mol_wt = rdMolDescriptors.CalcExactMolWt(product)

        return {
            "has_double_bond": has_double_bond,
            "product": product_smiles,
            "molecular_weight": round(mol_wt, 4)
        }
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "C=CC"
    print(f"result: {level_function(smiles)}")
