from rdkit import Chem

def level_function(mol):
    """
    将分子转为邻接矩阵。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        adj_matrix = Chem.GetAdjacencyMatrix(mol)
        return adj_matrix.tolist()
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CC[C@H](F)C(=O)O"
    print(f"邻接矩阵: {level_function(smiles)}")