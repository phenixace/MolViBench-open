from rdkit import Chem

def level_function(mol):
    """
    统计分子中卤素原子的数量。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        atoms = mol.GetAtoms()
        halogens = {"F", "Cl", "Br", "I", "At"}
        return sum(1 for atom in atoms if atom.GetSymbol() in halogens)
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "ClC(N)C[C@H](F)C(=O)O"
    print(f"卤素原子数量: {level_function(smiles)}")