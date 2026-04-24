from rdkit import Chem


# 常见毒性子结构 SMARTS
TOXIC_SMARTS = [
    "[$(a[N+](=O)[O-]),$(a[N](=O)=O)]",  # 硝基芳烃
    "[NH2]a",                              # 芳香胺
    "C1OC1",                               # 环氧化物
    "[CX3](=[OX1])[F,Cl,Br,I]",          # 酰卤
    "[N]=[N]",                             # 偶氮
    "OO",                                  # 过氧化物
    "[CX4]([F,Cl,Br,I])([F,Cl,Br,I])",   # 多卤代烃
]


def level_function(mols):
    """给定一组分子，过滤掉含有毒性子结构的分子。"""
    try:
        patterns = []
        for smarts in TOXIC_SMARTS:
            pat = Chem.MolFromSmarts(smarts)
            if pat:
                patterns.append(pat)

        results = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            is_toxic = False
            for pat in patterns:
                if mol.HasSubstructMatch(pat):
                    is_toxic = True
                    break
            if not is_toxic:
                results.append(Chem.MolToSmiles(mol))
        return results
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = ["CCO", "c1ccc([N+](=O)[O-])cc1", "c1ccccc1"]
    print(f"毒性过滤后: {level_function(smiles_list)}")
