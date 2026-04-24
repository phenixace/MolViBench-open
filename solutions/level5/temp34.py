from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


def level_function(mols):
    """给定一组分子，进行 Pareto 优化（目标：LogP ~2，TPSA < 120）。"""
    try:
        mol_data = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            logp = Descriptors.MolLogP(mol)
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            mol_data.append({
                'smiles': Chem.MolToSmiles(mol),
                'logp': round(logp, 2),
                'tpsa': round(tpsa, 2),
                'logp_dist': round(abs(logp - 2), 2)  # 离目标 LogP=2 的距离
            })

        # Pareto 前沿: 对两个目标 (logp_dist 越小越好, tpsa 越小越好)
        pareto_front = []
        for i, d in enumerate(mol_data):
            dominated = False
            for j, d2 in enumerate(mol_data):
                if i == j:
                    continue
                if (d2['logp_dist'] <= d['logp_dist'] and d2['tpsa'] <= d['tpsa'] and
                        (d2['logp_dist'] < d['logp_dist'] or d2['tpsa'] < d['tpsa'])):
                    dominated = True
                    break
            if not dominated and d['tpsa'] < 120:
                pareto_front.append(d)

        pareto_front.sort(key=lambda x: x['logp_dist'])
        return pareto_front
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "c1ccc(O)cc1",
                   "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "c1ccncc1"]
    result = level_function(smiles_list)
    if result:
        for r in result:
            print(f"  {r['smiles']}: LogP={r['logp']}, TPSA={r['tpsa']}")
