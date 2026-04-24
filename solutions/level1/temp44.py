from rdkit import Chem

def level_function(mol):
    """
    统计分子中氧原子的数量。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        atoms = mol.GetAtoms()
        oxygen_count = sum(1 for atom in atoms if atom.GetSymbol() == "O")
        return oxygen_count
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "PC(N)C[C@H](F)C(=O)O"
    print(f"氧原子数量: {level_function(smiles)}")