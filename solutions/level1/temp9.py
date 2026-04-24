from rdkit import Chem

def level_function(mol):
    """判断分子是否含有氨基 (-NH2)。"""
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return mol.HasSubstructMatch(Chem.MolFromSmarts("[NX3;H2]"))
    except Exception as e:
        print(e)
        return None
    
if __name__ == "__main__":
    smiles = "CC(N)C(=O)O"
    print(f"是否含有氨基: {level_function(smiles)}")
