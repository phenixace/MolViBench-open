from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


def level_function(mol):
    """给定分子 → 判断是否含醛基 → 若有 → 还原为醇 → 计算 LogP。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含醛基
        pattern = Chem.MolFromSmarts('[CX3H1](=O)')
        has_aldehyde = mol_obj.HasSubstructMatch(pattern)

        if not has_aldehyde:
            return None

        # Step 2: 还原醛基为醇
        rxn = AllChem.ReactionFromSmarts('[C:1](=O)[H]>>[C:1]O')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        # Step 3: 计算 LogP
        logp = rdMolDescriptors.CalcCrippenDescriptors(product)[0]

        return {
            "has_aldehyde": has_aldehyde,
            "product": product_smiles,
            "logp": round(logp, 4)
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "CC=O"
    print(f"result: {level_function(smiles)}")
