from rdkit import Chem


def level_function(filename):
    """从 SDF 文件中读取分子。"""
    try:
        supplier = Chem.SDMolSupplier(filename)
        smiles_list = []
        for mol in supplier:
            if mol is not None:
                smiles_list.append(Chem.MolToSmiles(mol))
        return smiles_list
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    # First create a test SDF file
    mol = Chem.MolFromSmiles("CCO")
    writer = Chem.SDWriter("test_input.sdf")
    writer.write(mol)
    writer.close()
    result = level_function("test_input.sdf")
    print(f"从 SDF 读取的分子: {result}")
