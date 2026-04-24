from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


def level_function(mol):
    """给定分子 → 判断是否含芳环 → 若有 → 甲基化 → 计算 TPSA。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含芳环
        num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol_obj)
        has_aromatic_ring = num_aromatic_rings > 0

        if not has_aromatic_ring:
            return None

        # Step 2: 甲基化 (芳环上添加甲基)
        rxn = AllChem.ReactionFromSmarts('[cH:1]>>[c:1]C')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        # Step 3: 计算 TPSA
        tpsa = rdMolDescriptors.CalcTPSA(product)

        return {
            "has_aromatic_ring": has_aromatic_ring,
            "product": product_smiles,
            "tpsa": round(tpsa, 4)
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccccc1O"
    print(f"result: {level_function(smiles)}")
