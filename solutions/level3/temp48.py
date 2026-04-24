from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


# REACH 相关环保规则 (简化版)
# 主要检查: PBT (持久性、生物蓄积性、毒性)
PBT_TOXIC_SMARTS = [
    "[$(a[N+](=O)[O-])]",  # 硝基芳烃
    "c1cc([Cl])c([Cl])c([Cl])c1",  # 多氯联苯类
    "[Sn]",  # 有机锡
    "[Pb]",  # 铅
    "[Hg]",  # 汞
    "[Cd]",  # 镉
]


def level_function(mol):
    """给定分子，预测是否符合 REACH 环保标准。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        issues = []

        # 1. 检查 LogP (高 LogP 可能意味着高生物蓄积性)
        logp = Descriptors.MolLogP(mol_obj)
        if logp > 4.5:
            issues.append(f"LogP={round(logp, 2)} > 4.5, 可能具有生物蓄积性")

        # 2. 检查分子量 (高分子量物质可能不易降解)
        mw = Descriptors.MolWt(mol_obj)

        # 3. 检查毒性子结构
        for smarts in PBT_TOXIC_SMARTS:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol_obj.HasSubstructMatch(pattern):
                issues.append(f"含有 PBT 相关毒性子结构: {smarts}")

        # 4. 检查卤素含量 (多卤化合物可能持久)
        halogen_count = sum(1 for atom in mol_obj.GetAtoms()
                           if atom.GetAtomicNum() in [9, 17, 35, 53])
        if halogen_count >= 3:
            issues.append(f"含有 {halogen_count} 个卤素原子, 可能具有持久性")

        return {
            "compliant": len(issues) == 0,
            "issues": issues,
            "logp": round(logp, 2),
            "mw": round(mw, 2)
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1cc(Cl)c(Cl)c(Cl)c1"
    print(f"REACH 评估: {level_function(smiles)}")
