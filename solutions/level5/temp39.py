from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors, Lipinski


def level_function(mols):
    """给定一组分子，找出相似度 <0.4 且满足 Lipinski 的分子（增加化学空间覆盖）。"""
    try:
        mol_data = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Lipinski.NumHDonors(mol)
            hba = Lipinski.NumHAcceptors(mol)
            passes_lipinski = mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10
            mol_data.append({
                'smiles': Chem.MolToSmiles(mol),
                'fp': fp,
                'passes_lipinski': passes_lipinski
            })

        if len(mol_data) < 2:
            return None

        # 找 pairwise 相似度 < 0.4 的分子对
        results = []
        for i, d in enumerate(mol_data):
            if not d['passes_lipinski']:
                continue
            is_diverse = True
            for j, d2 in enumerate(mol_data):
                if i == j:
                    continue
                sim = DataStructs.TanimotoSimilarity(d['fp'], d2['fp'])
                if sim >= 0.4:
                    is_diverse = False
                    break

            if is_diverse:
                results.append(d['smiles'])

        # 如果没有完全不相似的, 选最不相似的
        if not results:
            for d in mol_data:
                if d['passes_lipinski']:
                    avg_sim = sum(
                        DataStructs.TanimotoSimilarity(d['fp'], d2['fp'])
                        for d2 in mol_data if d2['smiles'] != d['smiles']
                    ) / (len(mol_data) - 1)
                    if avg_sim < 0.4:
                        results.append(d['smiles'])

        return results
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "c1ccncc1", "CCCCCCCC"]
    result = level_function(smiles_list)
    print(f"低相似度 + Lipinski: {result}")
