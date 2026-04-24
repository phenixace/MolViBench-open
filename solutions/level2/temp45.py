from rdkit import Chem


def level_function(mol1, mol2):
    """比较两个分子的规范化 SMILES 是否一致。"""
    try:
        m1 = Chem.MolFromSmiles(mol1)
        m2 = Chem.MolFromSmiles(mol2)
        if m1 is None or m2 is None:
            return False
        canonical1 = Chem.MolToSmiles(m1)
        canonical2 = Chem.MolToSmiles(m2)
        return canonical1 == canonical2
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    result = level_function("CCO", "OCC")
    print(f"两个分子的规范化 SMILES 是否一致: {result}")
    result2 = level_function("CCO", "CCCO")
    print(f"CCO 与 CCCO 是否一致: {result2}")
