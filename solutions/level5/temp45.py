from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors


def level_function(scaffold, n_modifications=10):
    """给定一个 scaffold → 生成多种侧链修饰 → 保留 QED > 0.7 → 比较相似度。"""
    try:
        scaffold_mol = Chem.MolFromSmiles(scaffold)
        if scaffold_mol is None:
            return None

        scaffold_fp = AllChem.GetMorganFingerprintAsBitVect(scaffold_mol, 2, nBits=2048)

        # Step 1: 生成侧链修饰
        rxns = [
            ('[cH:1]>>[c:1]C', '甲基'),
            ('[cH:1]>>[c:1]CC', '乙基'),
            ('[cH:1]>>[c:1]O', '羟基'),
            ('[cH:1]>>[c:1]OC', '甲氧基'),
            ('[cH:1]>>[c:1]F', '氟'),
            ('[cH:1]>>[c:1]Cl', '氯'),
            ('[cH:1]>>[c:1]N', '氨基'),
            ('[cH:1]>>[c:1]NC(=O)C', '酰胺'),
            ('[cH:1]>>[c:1]C(F)(F)F', '三氟甲基'),
            ('[cH:1]>>[c:1]C#N', '氰基'),
        ]

        derivatives = []
        seen = set()
        for rxn_smarts, name in rxns:
            rxn = AllChem.ReactionFromSmarts(rxn_smarts)
            products = rxn.RunReactants((scaffold_mol,))
            for product_set in products:
                for product in product_set:
                    try:
                        Chem.SanitizeMol(product)
                        smi = Chem.MolToSmiles(product)
                        if smi in seen:
                            continue
                        seen.add(smi)

                        # Step 2: 筛选 QED > 0.7
                        qed = Descriptors.qed(product)
                        if qed > 0.7:
                            fp = AllChem.GetMorganFingerprintAsBitVect(product, 2, nBits=2048)
                            sim = DataStructs.TanimotoSimilarity(scaffold_fp, fp)
                            derivatives.append({
                                'smiles': smi,
                                'qed': round(qed, 4),
                                'similarity_to_scaffold': round(sim, 4),
                                'modification': name
                            })
                    except Exception:
                        continue

        derivatives.sort(key=lambda x: x['qed'], reverse=True)
        return derivatives[:n_modifications]
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    scaffold = "c1ccccc1"
    result = level_function(scaffold)
    if result:
        for r in result[:5]:
            print(f"  {r['smiles']}: QED={r['qed']}, sim={r['similarity_to_scaffold']}")
