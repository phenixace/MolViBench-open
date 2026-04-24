from rdkit import Chem

def level_function(mol):
    """
    判断分子是否含有硫原子。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        atoms = mol.GetAtoms()
        for atom in atoms:
            if atom.GetSymbol() == "S":
                return True
        return False
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "SC(N)C[C@H](F)C(=O)O"
    print(f"是否含有硫原子: {level_function(smiles)}")