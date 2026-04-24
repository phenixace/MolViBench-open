from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):
    """给定分子 → 判断是否含炔基 → 若有 → 氢化为烯 → 计算分子量。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含炔基
        pattern = Chem.MolFromSmarts('[C:1]#[C:2]')
        has_alkyne = mol_obj.HasSubstructMatch(pattern)

        if not has_alkyne:
            return None

        # Step 2: 部分氢化（将 C≡C 变为 C=C）
        rxn = AllChem.ReactionFromSmarts('[C:1]#[C:2]>>[C:1]=[C:2]')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        # Step 3: 计算分子量
        mol_wt = rdMolDescriptors.CalcExactMolWt(product)

        return {
            "has_alkyne": has_alkyne,
            "product": product_smiles,
            "molecular_weight": round(mol_wt, 4)
        }
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CC#CC"
    print(f"result: {level_function(smiles)}")
