from rdkit import Chem


def level_function(mol, filename="output.mol"):
    """将分子导出为 MOL 文件。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        Chem.MolToMolFile(mol_obj, filename)
        return filename
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "CCO"
    result = level_function(smiles, "output.mol")
    print(f"分子已导出为 MOL 文件: {result}")
