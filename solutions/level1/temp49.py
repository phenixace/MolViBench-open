from rdkit import Chem

def level_function(mol):
    """
    统计分子中氢原子的数量。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)  # 添加显式氢原子
        atoms = mol.GetAtoms()
        atom_count = sum(1 for atom in atoms if atom.GetSymbol() == "H")
        
        return atom_count
        
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "SC(N)C[C@H](F)C(=O)O"
    print(f"氢原子数量: {level_function(smiles)}")
