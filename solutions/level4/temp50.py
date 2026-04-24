from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdMolDescriptors
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions


def level_function(mol):
    """给定分子 → 判断是否含苯环 → 若有 → 生成所有对映体 → 再计算相似度矩阵。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含苯环
        pattern = Chem.MolFromSmarts('c1ccccc1')
        has_benzene = mol_obj.HasSubstructMatch(pattern)

        if not has_benzene:
            return None

        # Step 2: 生成所有立体异构体
        opts = StereoEnumerationOptions(unique=True)
        isomers = list(EnumerateStereoisomers(mol_obj, options=opts))

        if len(isomers) < 1:
            return None

        isomer_smiles = []
        fps = []
        for iso in isomers:
            smi = Chem.MolToSmiles(iso)
            if smi not in isomer_smiles:
                isomer_smiles.append(smi)
                fp = AllChem.GetMorganFingerprintAsBitVect(iso, 2, nBits=2048)
                fps.append(fp)

        # Step 3: 计算相似度矩阵
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

        return {
            "has_benzene": has_benzene,
            "isomers": isomer_smiles,
            "similarity_matrix": matrix
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccc(CC(F)Cl)cc1"
    result = level_function(smiles)
    if result:
        print(f"异构体数量: {len(result['isomers'])}")
        for smi in result['isomers']:
            print(f"  {smi}")
