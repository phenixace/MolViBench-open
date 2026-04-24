from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol):
    """
    给定苯，生成所有一取代甲基衍生物。
    """
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        # 使用反应 SMARTS: 芳香碳上的氢替换为甲基
        rxn = AllChem.ReactionFromSmarts('[cH:1]>>[c:1]C')
        products = rxn.RunReactants((mol_obj,))
        # 收集所有唯一产物
        unique_smiles = set()
        for product_tuple in products:
            for product in product_tuple:
                Chem.SanitizeMol(product)
                smi = Chem.MolToSmiles(product)
                unique_smiles.add(smi)
        return sorted(list(unique_smiles))
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "c1ccccc1"  # 苯
    result = level_function(smiles)
    print(f"一取代甲基衍生物: {result}")
