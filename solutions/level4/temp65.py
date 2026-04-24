from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


def level_function(seed_smiles):
    """从种子分子出发 → 第一轮：生成所有单取代甲基衍生物 → 第二轮：对上轮最优分子生成所有单取代羟基衍生物 → 第三轮：对上轮最优分子生成所有单取代卤素衍生物 → 最终选 QED 最高的分子。"""
    try:
        seed = Chem.MolFromSmiles(seed_smiles)
        if seed is None:
            return None

        def generate_derivatives(mol, rxn_smarts):
            rxn = AllChem.ReactionFromSmarts(rxn_smarts)
            products = rxn.RunReactants((mol,))
            derivs = set()
            for prod_set in products:
                for prod in prod_set:
                    try:
                        Chem.SanitizeMol(prod)
                        derivs.add(Chem.MolToSmiles(prod))
                    except Exception:
                        pass
            return derivs

        def best_by_qed(smiles_set):
            best_smi = None
            best_qed = -1
            for smi in smiles_set:
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    q = Descriptors.qed(mol)
                    if q > best_qed:
                        best_qed = q
                        best_smi = smi
            return best_smi, best_qed

        rounds = []

        # Round 1: methyl derivatives
        methyl_derivs = generate_derivatives(seed, '[cH:1]>>[c:1]C')
        if not methyl_derivs:
            methyl_derivs = generate_derivatives(seed, '[CH:1]>>[C:1](C)')
        best1_smi, best1_qed = best_by_qed(methyl_derivs)
        rounds.append({"round": 1, "type": "methyl", "num_derivs": len(methyl_derivs),
                       "best": best1_smi, "QED": round(best1_qed, 4) if best1_qed > 0 else None})

        # Round 2: hydroxyl derivatives of best from round 1
        if best1_smi:
            mol1 = Chem.MolFromSmiles(best1_smi)
            oh_derivs = generate_derivatives(mol1, '[cH:1]>>[c:1]O')
            if not oh_derivs:
                oh_derivs = generate_derivatives(mol1, '[CH:1]>>[C:1](O)')
            best2_smi, best2_qed = best_by_qed(oh_derivs)
            rounds.append({"round": 2, "type": "OH", "num_derivs": len(oh_derivs),
                          "best": best2_smi, "QED": round(best2_qed, 4) if best2_qed > 0 else None})
        else:
            best2_smi = best1_smi
            best2_qed = best1_qed

        # Round 3: halogen derivatives of best from round 2
        if best2_smi:
            mol2 = Chem.MolFromSmiles(best2_smi)
            halogen_derivs = set()
            for hal in ['F', 'Cl', 'Br']:
                derivs = generate_derivatives(mol2, f'[cH:1]>>[c:1]{hal}')
                halogen_derivs.update(derivs)
            best3_smi, best3_qed = best_by_qed(halogen_derivs)
            rounds.append({"round": 3, "type": "halogen", "num_derivs": len(halogen_derivs),
                          "best": best3_smi, "QED": round(best3_qed, 4) if best3_qed > 0 else None})
        else:
            best3_smi = best2_smi
            best3_qed = best2_qed

        # Final best
        all_candidates = {seed_smiles}
        if best1_smi: all_candidates.add(best1_smi)
        if best2_smi: all_candidates.add(best2_smi)
        if best3_smi: all_candidates.add(best3_smi)
        final_smi, final_qed = best_by_qed(all_candidates)

        return {
            "seed": seed_smiles,
            "rounds": rounds,
            "final_best": final_smi,
            "final_QED": round(final_qed, 4) if final_qed > 0 else None
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    result = level_function("c1ccccc1")
    if result:
        print(f"Final: {result['final_best']}, QED: {result['final_QED']}")
        for r in result['rounds']:
            print(f"  Round {r['round']}: {r['type']}, {r['num_derivs']} derivs, best QED={r['QED']}")
