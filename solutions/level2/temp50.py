from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs


def level_function(mols):
    """计算一组分子的 pairwise 相似度矩阵。"""
    try:
        fps = []
        valid_smiles = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fps.append(fp)
            valid_smiles.append(smi)
        if len(fps) < 2:
            return None

        n = len(fps)
        matrix = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i == j:
                    matrix[i][j] = 1.0
                elif j > i:
                    sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
                    matrix[i][j] = round(sim, 4)
                    matrix[j][i] = round(sim, 4)

        return {"smiles": valid_smiles, "similarity_matrix": matrix}
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "CCCC"]
    result = level_function(smiles_list)
    if result:
        print("SMILES:", result["smiles"])
        for row in result["similarity_matrix"]:
            print([f"{x:.4f}" for x in row])
