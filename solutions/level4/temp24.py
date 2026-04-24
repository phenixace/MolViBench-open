from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):
    """给定分子 → 判断是否含氨基 → 若有 → 烷基化 → 计算 LogP。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含氨基
        pattern = Chem.MolFromSmarts('[NX3;H2,H1]')
        has_amino = mol_obj.HasSubstructMatch(pattern)

        if not has_amino:
            return None

        # Step 2: N-甲基化
        rxn = AllChem.ReactionFromSmarts('[NH2:1]>>[NH:1]C')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        # Step 3: 计算 LogP
        logp = rdMolDescriptors.CalcCrippenDescriptors(product)[0]

        return {
            "has_amino": has_amino,
            "product": product_smiles,
            "logp": round(logp, 4)
        }
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "c1ccc(N)cc1"
    print(f"result: {level_function(smiles)}")
