from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors


def level_function(mols):
    """给定一组分子，构建"相似度–药物性"热图。"""
    try:
        mol_data = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            qed = Descriptors.qed(mol)
            mol_data.append({
                'smiles': Chem.MolToSmiles(mol),
                'fp': fp,
                'qed': round(qed, 4)
            })

        n = len(mol_data)
        if n < 2:
            return None

        # 构建相似度矩阵
        sim_matrix = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i == j:
                    sim_matrix[i][j] = 1.0
                elif j > i:
                    sim = DataStructs.TanimotoSimilarity(
                        mol_data[i]['fp'], mol_data[j]['fp']
                    )
                    sim_matrix[i][j] = round(sim, 4)
                    sim_matrix[j][i] = round(sim, 4)

        return {
            'smiles': [d['smiles'] for d in mol_data],
            'qed_values': [d['qed'] for d in mol_data],
            'similarity_matrix': sim_matrix
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "c1ccncc1"]
    result = level_function(smiles_list)
    if result:
        print(f"SMILES: {result['smiles']}")
        print(f"QED: {result['qed_values']}")
        print("相似度矩阵:")
        for row in result['similarity_matrix']:
            print(f"  {row}")
