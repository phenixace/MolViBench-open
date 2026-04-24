from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors


def level_function(mol):
    """从一个起始分子 → 生成所有可能的衍生物 → 多目标优化（QED, LogP, TPSA）→ 选出最佳 lead。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        orig_smi = Chem.MolToSmiles(mol_obj)

        # Step 1: 生成衍生物
        rxns = [
            '[cH:1]>>[c:1]C',
            '[cH:1]>>[c:1]O',
            '[cH:1]>>[c:1]N',
            '[cH:1]>>[c:1]F',
            '[cH:1]>>[c:1]Cl',
            '[cH:1]>>[c:1]OC',
            '[cH:1]>>[c:1]C(F)(F)F',
            '[cH:1]>>[c:1]C#N',
            '[cH:1]>>[c:1]CC',
            '[cH:1]>>[c:1]NC(=O)C',
        ]

        candidates = set()
        candidates.add(orig_smi)

        for rxn_smarts in rxns:
            rxn = AllChem.ReactionFromSmarts(rxn_smarts)
            products = rxn.RunReactants((mol_obj,))
            for product_set in products:
                for product in product_set:
                    try:
                        Chem.SanitizeMol(product)
                        smi = Chem.MolToSmiles(product)
                        candidates.add(smi)

                        # 二次修饰
                        products2 = rxn.RunReactants((product,))
                        for ps2 in products2:
                            for p2 in ps2:
                                try:
                                    Chem.SanitizeMol(p2)
                                    smi2 = Chem.MolToSmiles(p2)
                                    candidates.add(smi2)
                                except Exception:
                                    continue
                    except Exception:
                        continue

        # Step 2: 多目标优化 (QED 高, LogP ∈ [1,3], TPSA < 120)
        scored = []
        for smi in candidates:
            m = Chem.MolFromSmiles(smi)
            if m is None:
                continue
            qed = Descriptors.qed(m)
            logp = Descriptors.MolLogP(m)
            tpsa = rdMolDescriptors.CalcTPSA(m)

            # 综合评分
            logp_score = 1.0 - min(abs(logp - 2.0), 3.0) / 3.0  # LogP 越接近 2 越好
            tpsa_score = 1.0 if tpsa < 120 else max(0, 1.0 - (tpsa - 120) / 60.0)
            total_score = qed * 0.4 + logp_score * 0.3 + tpsa_score * 0.3

            scored.append({
                'smiles': smi,
                'qed': round(qed, 4),
                'logp': round(logp, 2),
                'tpsa': round(tpsa, 2),
                'total_score': round(total_score, 4)
            })

        # Step 3: 选出最佳 lead
        scored.sort(key=lambda x: x['total_score'], reverse=True)

        return {
            'num_candidates': len(candidates),
            'best_lead': scored[0] if scored else None,
            'top5': scored[:5]
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccccc1"
    result = level_function(smiles)
    if result:
        print(f"候选数量: {result['num_candidates']}")
        print(f"最佳 Lead: {result['best_lead']}")
        print("Top 5:")
        for r in result['top5']:
            print(f"  {r['smiles']}: score={r['total_score']}, QED={r['qed']}, "
                  f"LogP={r['logp']}, TPSA={r['tpsa']}")
