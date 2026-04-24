from rdkit import Chem

def level_function(mol):
    """
    检查分子 SMILES 是否有效。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return False
        return True
    except Exception as e:
        print(e)
        return False

if __name__ == "__main__":
    smiles = "CC[C@H](F)C(=O)O"
    print(f"分子 SMILES 是否有效: {level_function(smiles)}")