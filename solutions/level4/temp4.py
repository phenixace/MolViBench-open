from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, RWMol

def level_function(mol):
    """给定分子 → 判断是否含卤素 → 若有 → 替换为氢 → 计算 LogP。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含卤素
        pattern = Chem.MolFromSmarts('[F,Cl,Br,I]')
        has_halogen = mol_obj.HasSubstructMatch(pattern)

        if not has_halogen:
            return None

        # Step 2: 替换卤素为氢
        rxn = AllChem.ReactionFromSmarts('[C:1][F,Cl,Br,I]>>[C:1][H]')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        # Step 3: 计算 LogP
        logp = rdMolDescriptors.CalcCrippenDescriptors(product)[0]

        return {
            "has_halogen": has_halogen,
            "product": product_smiles,
            "logp": round(logp, 4)
        }
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "c1ccc(Cl)cc1"
    print(f"result: {level_function(smiles)}")
