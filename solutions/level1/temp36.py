from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def level_function(mol1, mol2, consider_stereochemistry=True):
    """
    计算两个分子是否为同分异构体。
    
    参数:
        mol1, mol2: SMILES 字符串
        consider_stereochemistry: 是否考虑立体异构 (默认 True)
    返回:
        True / False / None
    """
    try:
        mol1 = Chem.MolFromSmiles(mol1)
        mol2 = Chem.MolFromSmiles(mol2)
        if mol1 is None or mol2 is None:
            return None

        # 1. 判断分子式是否相同
        formula1 = rdMolDescriptors.CalcMolFormula(mol1)
        formula2 = rdMolDescriptors.CalcMolFormula(mol2)
        if formula1 != formula2:
            return False  # 不是同分异构体

        # 2. 获取规范化 SMILES
        smiles1 = Chem.MolToSmiles(mol1, isomericSmiles=consider_stereochemistry, canonical=True)
        smiles2 = Chem.MolToSmiles(mol2, isomericSmiles=consider_stereochemistry, canonical=True)

        if smiles1 != smiles2:
            return True   # 分子式相同但结构不同 → 同分异构体
        else:
            return False  # 完全相同

    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles1 = "CC[C@H](F)C(=O)O"   # 有手性
    smiles2 = "CC[C@@H](F)C(=O)O"  # 对映体
    smiles3 = "CCC(F)C(=O)O"       # 去掉立体信息
    
    print("考虑立体化学：")
    print(f"{smiles1} vs {smiles2} → {level_function(smiles1, smiles2, consider_stereochemistry=True)}")
    print(f"{smiles1} vs {smiles3} → {level_function(smiles1, smiles3, consider_stereochemistry=True)}")

    print("\n不考虑立体化学：")
    print(f"{smiles1} vs {smiles2} → {level_function(smiles1, smiles2, consider_stereochemistry=False)}")
    print(f"{smiles1} vs {smiles3} → {level_function(smiles1, smiles3, consider_stereochemistry=False)}")
