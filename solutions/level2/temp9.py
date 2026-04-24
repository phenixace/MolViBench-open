from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol, substituent="C"):
    """
    给定苯，生成对位取代产物。
    """
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        # 对位取代：在苯环的1,4位添加取代基
        # 使用反应 SMARTS 进行对位取代
        rxn_smarts = f'[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1>>[c:1]({substituent})1[cH:2][cH:3][c:4]({substituent})[cH:5][cH:6]1'
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
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
    smiles = "c1ccccc1"  # 苯
    result = level_function(smiles, "C")
    print(f"对位取代产物: {result}")
