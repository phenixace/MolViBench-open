from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors


def level_function(mol):
    """给定一个候选分子，生成更极性的衍生物以提高溶解度。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        orig_tpsa = rdMolDescriptors.CalcTPSA(mol_obj)
        orig_smi = Chem.MolToSmiles(mol_obj)

        # 策略: 添加极性基团 (-OH, -NH2, =O)
        polar_rxns = [
            ('[cH:1]>>[c:1]O', '添加羟基'),
            ('[cH:1]>>[c:1]N', '添加氨基'),
            ('[CH3:1]>>[CH2:1]O', '烷基氧化'),
        ]

        derivatives = []
        for smarts, desc in polar_rxns:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((mol_obj,))
            for product_set in products:
                for product in product_set:
                    try:
                        Chem.SanitizeMol(product)
                        smi = Chem.MolToSmiles(product)
                        if smi == orig_smi or smi in [d['smiles'] for d in derivatives]:
                            continue
                        new_tpsa = rdMolDescriptors.CalcTPSA(product)
                        if new_tpsa > orig_tpsa:
                            derivatives.append({
                                'smiles': smi,
                                'tpsa': round(new_tpsa, 2),
                                'modification': desc
                            })
                    except Exception:
                        continue

        derivatives.sort(key=lambda x: x['tpsa'], reverse=True)
        return derivatives[:10]
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccccc1CC"
    result = level_function(smiles)
    print(f"更极性的衍生物: {result}")
