from rdkit import Chem


# 常见毒性子结构 SMARTS (简化版)
TOXIC_SUBSTRUCTURES = {
    "硝基芳烃": "[$(a[N+](=O)[O-]),$(a[N](=O)=O)]",
    "芳香胺": "[NH2]a",
    "醛基": "[CH1](=O)",
    "环氧化物": "C1OC1",
    "酰卤": "[CX3](=[OX1])[F,Cl,Br,I]",
    "异氰酸酯": "[N]=[C]=[O]",
    "偶氮化合物": "[N]=[N]",
    "亚硝基": "[N]=O",
    "磺酸酯": "S(=O)(=O)O[C,c]",
    "磷酸酯": "P(=O)(O)(O)O",
    "过氧化物": "OO",
    "迈克尔受体": "[CH2]=[CH][C,S,N](=O)",
    "卤代烃 (多卤)": "[CX4]([F,Cl,Br,I])([F,Cl,Br,I])",
}


def level_function(mol):
    """给定分子，预测可能的毒性子结构。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        found_alerts = []
        for name, smarts in TOXIC_SUBSTRUCTURES.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                continue
            if mol_obj.HasSubstructMatch(pattern):
                found_alerts.append(name)

        return found_alerts if found_alerts else []
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccc([N+](=O)[O-])cc1N"  # 对硝基苯胺
    print(f"毒性子结构: {level_function(smiles)}")
