from rdkit import Chem

def level_function(mol):
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        adj_matrix = Chem.GetAdjacencyMatrix(mol)
        return adj_matrix.tolist()
    except Exception as e:
        print(e)
        return None
