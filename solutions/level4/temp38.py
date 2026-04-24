from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


def level_function(mol):
    """给定分子 → 判断是否含胺基 → 若有 → 磺酰化 → 计算 TPSA。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含胺基
        pattern = Chem.MolFromSmarts('[NX3;H2,H1]')
        has_amine = mol_obj.HasSubstructMatch(pattern)

        if not has_amine:
            return None

        # Step 2: 磺酰化 (R-NH2 → R-NHS(=O)(=O)C)
        rxn = AllChem.ReactionFromSmarts('[NH2:1]>>[NH:1]S(=O)(=O)C')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            # 尝试 -NHR
            rxn = AllChem.ReactionFromSmarts('[NH1:1]>>[N:1]S(=O)(=O)C')
            products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        # Step 3: 计算 TPSA
        tpsa = rdMolDescriptors.CalcTPSA(product)

        return {
            "has_amine": has_amine,
            "product": product_smiles,
            "tpsa": round(tpsa, 4)
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "CCN"
    print(f"result: {level_function(smiles)}")
