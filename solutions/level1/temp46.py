from rdkit import Chem

def level_function(mol):
    """
    统计分子中碳原子的数量。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        atoms = mol.GetAtoms()
        return sum(1 for atom in atoms if atom.GetSymbol() == "C")
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "PC(N)C[C@H](F)C(=O)O"
    print(f"碳原子数量: {level_function(smiles)}")