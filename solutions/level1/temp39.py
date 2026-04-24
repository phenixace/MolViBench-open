from rdkit import Chem

def level_function(mol):
    """
    判断分子是否含有芳香氮。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        pattern = Chem.MolFromSmarts("n")
        return mol.HasSubstructMatch(pattern)
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "C(N)C[C@H](F)C(=O)O"
    print(f"是否含有芳香氮: {level_function(smiles)}")