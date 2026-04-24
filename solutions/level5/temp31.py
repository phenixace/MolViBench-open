from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors


def level_function(mols):
    """给定一组分子，计算它们的 "相似度 vs QED" 二维图。"""
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

        if len(mol_data) < 2:
            return None

        # 计算每个分子与其他分子的平均相似度
        results = []
        for i, d in enumerate(mol_data):
            sims = []
            for j, d2 in enumerate(mol_data):
                if i != j:
                    sim = DataStructs.TanimotoSimilarity(d['fp'], d2['fp'])
                    sims.append(sim)
            avg_sim = sum(sims) / len(sims) if sims else 0
            results.append({
                'smiles': d['smiles'],
                'avg_similarity': round(avg_sim, 4),
                'qed': d['qed']
            })

        return results
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", "c1ccncc1"]
    result = level_function(smiles_list)
    if result:
        for r in result:
            print(f"  {r['smiles']}: sim={r['avg_similarity']}, QED={r['qed']}")
