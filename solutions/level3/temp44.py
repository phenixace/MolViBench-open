from rdkit import Chem
from rdkit.Chem import AllChem


# 常见代谢反应 SMARTS
METABOLIC_REACTIONS = [
    # 氧化反应
    ("[CH3:1]>>[CH2:1]O", "烷基氧化 (C-H → C-OH)"),
    ("[cH:1]>>[c:1]O", "芳环氧化"),
    ("[NH2:1]>>[NH:1]O", "N-氧化"),
    ("[SX2:1]>>[S:1](=O)", "S-氧化"),
    # 脱甲基反应
    ("[O:1]C>>[O:1]", "O-脱甲基"),
    ("[N:1](C)>>[NH:1]", "N-脱甲基"),
    # 水解
    ("[C:1](=O)[O:2][C:3]>>[C:1](=O)[OH:2].[OH][C:3]", "酯水解"),
    # 还原
    ("[C:1](=O)[C:2]>>[C:1](O)[C:2]", "酮还原"),
]


def level_function(mol):
    """给定分子，预测可能的代谢产物。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        metabolites = []
        for smarts, name in METABOLIC_REACTIONS:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((mol_obj,))
            for product_set in products:
                for product in product_set:
                    try:
                        Chem.SanitizeMol(product)
                        smi = Chem.MolToSmiles(product)
                        if smi not in [m["smiles"] for m in metabolites]:
                            metabolites.append({"smiles": smi, "reaction": name})
                    except Exception:
                        continue

        return metabolites if metabolites else None
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "CC(=O)Oc1ccccc1OC(C)=O"  # 阿司匹林
    result = level_function(smiles)
    if result:
        for m in result:
            print(f"  {m['reaction']}: {m['smiles']}")
