from rdkit import Chem

def level_function(mol):
    """
    获取分子的原子数。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return mol.GetNumAtoms()
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CCO"
    print(f"原子数: {level_function(smiles)}")