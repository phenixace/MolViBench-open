from rdkit import Chem


def level_function(mol):
    """给定分子，生成 SMILES 的所有规范化形式。"""
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        smiles_set = set()
        # Generate multiple random SMILES representations
        for _ in range(100):
            smi = Chem.MolToRandomSmilesVect(mol, 1)[0]
            smiles_set.add(smi)
        # Also add the canonical form
        smiles_set.add(Chem.MolToSmiles(mol))
        return list(smiles_set)
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    result = level_function("c1ccccc1CCO")
    print(f"SMILES 的多种表示 ({len(result)} 种): {result[:5]}...")
