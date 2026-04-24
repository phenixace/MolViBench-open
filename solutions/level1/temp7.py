from rdkit import Chem

def level_function(mol):
    """判断分子是否含有羟基 (-OH)。"""
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return mol.HasSubstructMatch(Chem.MolFromSmarts("[OX2H]"))
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CC(=O)O"
    print(f"是否含有羟基: {level_function(smiles)}")