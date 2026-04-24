from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):
    """给定分子 → 判断是否含芳香氮 → 若有 → 甲基化 → 计算 LogP。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含芳香氮
        pattern = Chem.MolFromSmarts('[n]')
        has_aromatic_n = mol_obj.HasSubstructMatch(pattern)

        if not has_aromatic_n:
            return None

        # Step 2: N-甲基化
        rxn = AllChem.ReactionFromSmarts('[n:1]>>[n+:1]C')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        # Step 3: 计算 LogP
        logp = rdMolDescriptors.CalcCrippenDescriptors(product)[0]

        return {
            "has_aromatic_n": has_aromatic_n,
            "product": product_smiles,
            "logp": round(logp, 4)
        }
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "c1ccncc1"
    print(f"result: {level_function(smiles)}")
