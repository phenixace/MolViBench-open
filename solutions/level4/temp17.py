from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):
    """给定分子 → 判断是否含酯基 → 若有 → 水解为羧酸 + 醇 → 计算分子式。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含酯基
        pattern = Chem.MolFromSmarts('C(=O)O[#6]')
        has_ester = mol_obj.HasSubstructMatch(pattern)

        if not has_ester:
            return None

        # Step 2: 水解酯基（生成羧酸 + 醇）
        rxn = AllChem.ReactionFromSmarts('[C:1](=O)[O:2][C:3]>>[C:1](=O)[OH].[OH][C:3]')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product1 = products[0][0]
        product2 = products[0][1]
        Chem.SanitizeMol(product1)
        Chem.SanitizeMol(product2)
        product1_smiles = Chem.MolToSmiles(product1)
        product2_smiles = Chem.MolToSmiles(product2)

        # Step 3: 计算两个产物的分子式
        formula1 = rdMolDescriptors.CalcMolFormula(product1)
        formula2 = rdMolDescriptors.CalcMolFormula(product2)

        return {
            "has_ester": has_ester,
            "product_acid": product1_smiles,
            "product_alcohol": product2_smiles,
            "formula_acid": formula1,
            "formula_alcohol": formula2
        }
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CC(=O)OCC"
    print(f"result: {level_function(smiles)}")
