from rdkit import Chem

def level_function(mol):
    """
    输出分子中所有原子符号。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return set([atom.GetSymbol() for atom in mol.GetAtoms()])
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CCOc1ccccc1"
    print(f"所有原子符号: {level_function(smiles)}")