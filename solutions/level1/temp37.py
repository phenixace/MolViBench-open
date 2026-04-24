from rdkit import Chem

def level_function(mol, substructure):
    """
    从分子中提取子结构匹配坐标。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        substructure = Chem.MolFromSmiles(substructure)
        if substructure is None:
            return None
        return mol.GetSubstructMatch(substructure)
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CC[C@H](F)C(=O)O"
    print(f"从分子中提取子结构匹配坐标。: {level_function(smiles, 'C(=O)O')}")