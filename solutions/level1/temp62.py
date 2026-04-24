from rdkit import Chem


def level_function(mol):
    """输出分子的 SMARTS 模式表示。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        return Chem.MolToSmarts(mol_obj)
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccccc1"
    print(f"SMARTS: {level_function(smiles)}")
