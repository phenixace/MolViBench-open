from rdkit import Chem
from rdkit.Chem import rdFMCS


def level_function(mol1, mol2):
    """计算两个分子的最大公共子结构（MCS）。"""
    try:
        m1 = Chem.MolFromSmiles(mol1)
        m2 = Chem.MolFromSmiles(mol2)
        if m1 is None or m2 is None:
            return None
        mcs_result = rdFMCS.FindMCS([m1, m2],
                                     atomCompare=rdFMCS.AtomCompare.CompareElements,
                                     bondCompare=rdFMCS.BondCompare.CompareOrder,
                                     timeout=10)
        return {
            "smartsString": mcs_result.smartsString,
            "numAtoms": mcs_result.numAtoms,
            "numBonds": mcs_result.numBonds,
            "canceled": mcs_result.canceled
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smi1 = "c1ccccc1CCO"
    smi2 = "c1ccccc1CCN"
    print(f"MCS: {level_function(smi1, smi2)}")
