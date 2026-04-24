from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol):
    """
    将甲苯中的甲基替换为羟基。
    """
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        # 使用反应 SMARTS: 将芳香碳上的甲基替换为羟基
        rxn = AllChem.ReactionFromSmarts('[c:1][CH3]>>[c:1]O')
        products = rxn.RunReactants((mol_obj,))
        unique_smiles = set()
        for product_tuple in products:
            for product in product_tuple:
                Chem.SanitizeMol(product)
                smi = Chem.MolToSmiles(product)
                unique_smiles.add(smi)
        if unique_smiles:
            return sorted(list(unique_smiles))[0]
        return None
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "Cc1ccccc1"  # 甲苯
    result = level_function(smiles)
    print(f"甲基替换为羟基: {result}")
