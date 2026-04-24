from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


def level_function(mol):
    """给定分子 → 判断是否含羟基 → 若有 → 卤代 → 再计算 QED。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含羟基
        pattern = Chem.MolFromSmarts('[OX2H]')
        has_hydroxyl = mol_obj.HasSubstructMatch(pattern)

        if not has_hydroxyl:
            return None

        # Step 2: 卤代 (替换 -OH 为 -Cl)
        rxn = AllChem.ReactionFromSmarts('[C:1][OH]>>[C:1]Cl')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        # Step 3: 计算 QED
        qed = Descriptors.qed(product)

        return {
            "has_hydroxyl": has_hydroxyl,
            "product": product_smiles,
            "qed": round(qed, 4)
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccc(O)cc1"
    print(f"result: {level_function(smiles)}")
