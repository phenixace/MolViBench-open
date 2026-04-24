from rdkit import Chem

# 常见官能团的 SMARTS 定义
FUNCTIONAL_GROUPS = {
    "羟基 (-OH)": "[OX2H]", 
    "氨基 (-NH2)": "[NX3;H2]", 
    "羧基 (-COOH)": "C(=O)[OX2H1]", 
    "醛基 (-CHO)": "[CX3H1](=O)[#6]", 
    "酮基 (C=O)": "[CX3](=O)[#6]", 
    "酯基 (-COOR)": "C(=O)O[#6]", 
    "卤素 (F/Cl/Br/I)": "[F,Cl,Br,I]", 
    "芳环 (Ar)": "a", 
    "膦基 (P)": "[PX4]", 
    "磺酰基 (-SO2-)": "S(=O)(=O)[#6]"
}

def level_function(mol):
    """
    获取分子的官能团列表。
    """
    try:
        mol = Chem.MolFromSmiles(mol) if isinstance(mol, str) else mol
        if mol is None:
            return None
        
        found_groups = []
        for name, smarts in FUNCTIONAL_GROUPS.items():
            patt = Chem.MolFromSmarts(smarts)
            if mol.HasSubstructMatch(patt):
                found_groups.append(name)
        
        return found_groups if found_groups else None

    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "PC(N)C[C@H](F)C(=O)O"  # 有 P, 氨基, F, 羧基
    print(f"官能团列表: {level_function(smiles)}")
