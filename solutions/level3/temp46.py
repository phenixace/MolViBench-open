from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


# 常见药效团 SMARTS 定义
PHARMACOPHORE_PATTERNS = {
    "氢键供体": "[#7H,#8H,#16H]",
    "氢键受体": "[#7,#8,#16;!H0;v2,v3,v4,v5]",
    "正电荷中心": "[+,NH3+,NH2+,NH+]",
    "负电荷中心": "[-,C(=O)[O-],S(=O)(=O)[O-]]",
    "芳香环": "a1aaaaa1",
    "疏水中心": "[CH2,CH3,c]",
    "卤素": "[F,Cl,Br,I]",
}


def level_function(mol):
    """给定分子，预测可能的药效团。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        pharmacophores = []
        for name, smarts in PHARMACOPHORE_PATTERNS.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                continue
            matches = mol_obj.GetSubstructMatches(pattern)
            if matches:
                pharmacophores.append({
                    "type": name,
                    "count": len(matches),
                    "atom_indices": [list(m) for m in matches]
                })

        return pharmacophores if pharmacophores else None
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "CC(=O)Oc1ccccc1C(=O)O"  # 阿司匹林
    result = level_function(smiles)
    if result:
        for p in result:
            print(f"  {p['type']}: {p['count']} 个")
