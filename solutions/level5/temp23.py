from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


def level_function(mol):
    """给定一个候选分子，优化其 LogP 使其落入 2~3。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        orig_logp = Descriptors.MolLogP(mol_obj)
        orig_smi = Chem.MolToSmiles(mol_obj)

        # 根据当前 LogP 选择策略
        if orig_logp < 2:
            # 需要增加 LogP: 添加疏水基团
            rxns = [
                ('[cH:1]>>[c:1]C', '添加甲基'),
                ('[cH:1]>>[c:1]F', '添加氟'),
                ('[cH:1]>>[c:1]Cl', '添加氯'),
            ]
        elif orig_logp > 3:
            # 需要降低 LogP: 添加极性基团
            rxns = [
                ('[cH:1]>>[c:1]O', '添加羟基'),
                ('[cH:1]>>[c:1]N', '添加氨基'),
                ('[CH3:1]>>[CH2:1]O', '烷基氧化'),
            ]
        else:
            return {"smiles": orig_smi, "logp": round(orig_logp, 2),
                    "message": "LogP 已在 2~3 范围内"}

        candidates = []
        for smarts, desc in rxns:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((mol_obj,))
            for product_set in products:
                for product in product_set:
                    try:
                        Chem.SanitizeMol(product)
                        smi = Chem.MolToSmiles(product)
                        if smi == orig_smi:
                            continue
                        logp = Descriptors.MolLogP(product)
                        if 2 <= logp <= 3:
                            candidates.append({
                                'smiles': smi,
                                'logp': round(logp, 2),
                                'modification': desc
                            })
                    except Exception:
                        continue

        # 选择最接近 2.5 的
        candidates.sort(key=lambda x: abs(x['logp'] - 2.5))
        return candidates[:5] if candidates else None
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccccc1"
    print(f"LogP 优化: {level_function(smiles)}")
