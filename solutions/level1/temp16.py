from rdkit import Chem

def level_function(mol):
    """
    获取分子的键数。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return mol.GetNumBonds()
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CCO"
    print(f"化学键数: {level_function(smiles)}")