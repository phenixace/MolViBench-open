from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


def level_function(mols):
    """给定一组分子，找出 QED 和 TPSA 最优的前 5 个分子。"""
    try:
        mol_data = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            qed = Descriptors.qed(mol)
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            mol_data.append({
                'smiles': Chem.MolToSmiles(mol),
                'qed': round(qed, 4),
                'tpsa': round(tpsa, 2)
            })

        # 排序: QED 越高越好, TPSA 越低越好 → 综合排名
        # 使用简单加权: score = QED - TPSA/200
        for d in mol_data:
            d['score'] = round(d['qed'] - d['tpsa'] / 200, 4)

        mol_data.sort(key=lambda x: x['score'], reverse=True)
        return mol_data[:5]
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O",
                   "c1ccncc1", "c1ccc(O)cc1", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"]
    result = level_function(smiles_list)
    if result:
        for r in result:
            print(f"  {r['smiles']}: QED={r['qed']}, TPSA={r['tpsa']}")
