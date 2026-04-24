from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors


# 简化药效团 SMARTS
PHARMACOPHORE_SMARTS = {
    "氢键供体": "[#7H,#8H]",
    "氢键受体": "[#7,#8]",
    "芳香环": "a1aaaaa1",
    "疏水": "[CH2,CH3]",
}


def level_function(pharmacophore_type="芳香环"):
    """从一个药效团 → 生成匹配分子 → 筛选 TPSA < 90 → 选择分子量最小的前 3 个。"""
    try:
        # Step 1: 基于药效团类型生成候选分子
        if pharmacophore_type not in PHARMACOPHORE_SMARTS:
            return None

        smarts = PHARMACOPHORE_SMARTS[pharmacophore_type]

        # 预定义一些含该药效团的分子骨架
        base_molecules = [
            "c1ccccc1", "c1ccncc1", "c1ccoc1", "c1ccsc1",
            "c1ccc(O)cc1", "c1ccc(N)cc1", "c1ccc(F)cc1",
            "c1ccc(C)cc1", "c1ccc(CC)cc1", "c1ccc2ccccc2c1",
        ]

        # 生成衍生物
        rxns = [
            '[cH:1]>>[c:1]C',
            '[cH:1]>>[c:1]O',
            '[cH:1]>>[c:1]F',
            '[cH:1]>>[c:1]N',
        ]

        candidates = set()
        for base_smi in base_molecules:
            base = Chem.MolFromSmiles(base_smi)
            if base is None:
                continue
            pattern = Chem.MolFromSmarts(smarts)
            if not base.HasSubstructMatch(pattern):
                continue
            candidates.add(base_smi)
            for rxn_smarts in rxns:
                rxn = AllChem.ReactionFromSmarts(rxn_smarts)
                products = rxn.RunReactants((base,))
                for ps in products:
                    for p in ps:
                        try:
                            Chem.SanitizeMol(p)
                            smi = Chem.MolToSmiles(p)
                            candidates.add(smi)
                        except Exception:
                            continue

        # Step 2: 筛选 TPSA < 90
        filtered = []
        for smi in candidates:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            if tpsa < 90:
                mw = Descriptors.MolWt(mol)
                filtered.append({
                    'smiles': smi,
                    'tpsa': round(tpsa, 2),
                    'mw': round(mw, 2)
                })

        # Step 3: 选择分子量最小的前 3 个
        filtered.sort(key=lambda x: x['mw'])
        return filtered[:3]
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    result = level_function("芳香环")
    if result:
        for r in result:
            print(f"  {r['smiles']}: MW={r['mw']}, TPSA={r['tpsa']}")
