from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


def level_function(mol):
    """给定一个候选分子，生成更大但保持 QED > 0.5 的衍生物。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        orig_mw = Descriptors.MolWt(mol_obj)
        orig_smi = Chem.MolToSmiles(mol_obj)

        # 策略: 添加基团
        growth_rxns = [
            ('[cH:1]>>[c:1]C', '添加甲基'),
            ('[cH:1]>>[c:1]OC', '添加甲氧基'),
            ('[cH:1]>>[c:1]F', '添加氟'),
            ('[cH:1]>>[c:1]Cl', '添加氯'),
            ('[cH:1]>>[c:1]CC', '添加乙基'),
            ('[NH2:1]>>[NH:1]C(=O)C', '酰化'),
        ]

        derivatives = []
        for smarts, desc in growth_rxns:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((mol_obj,))
            for product_set in products:
                for product in product_set:
                    try:
                        Chem.SanitizeMol(product)
                        smi = Chem.MolToSmiles(product)
                        if smi == orig_smi or smi in [d['smiles'] for d in derivatives]:
                            continue
                        new_mw = Descriptors.MolWt(product)
                        qed = Descriptors.qed(product)
                        if new_mw > orig_mw and qed > 0.5:
                            derivatives.append({
                                'smiles': smi,
                                'mw': round(new_mw, 2),
                                'qed': round(qed, 4),
                                'modification': desc
                            })
                    except Exception:
                        continue

        derivatives.sort(key=lambda x: x['qed'], reverse=True)
        return derivatives[:10]
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccc(N)cc1"
    result = level_function(smiles)
    print(f"更大的衍生物: {result}")
