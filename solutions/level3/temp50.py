from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


def level_function(mol):
    """给定分子，预测是否符合 Rule of Three（片段库）。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Rule of Three 条件:
        # MW ≤ 300
        # LogP ≤ 3
        # HBD ≤ 3
        # HBA ≤ 3
        # 旋转键 ≤ 3
        # TPSA ≤ 60

        mw = Descriptors.MolWt(mol_obj)
        logp = Descriptors.MolLogP(mol_obj)
        hbd = rdMolDescriptors.CalcNumHBD(mol_obj)
        hba = rdMolDescriptors.CalcNumHBA(mol_obj)
        rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol_obj)
        tpsa = rdMolDescriptors.CalcTPSA(mol_obj)

        passes = (mw <= 300 and logp <= 3 and hbd <= 3 and
                  hba <= 3 and rot_bonds <= 3 and tpsa <= 60)

        return {
            "passes_ro3": passes,
            "MW": round(mw, 2),
            "LogP": round(logp, 2),
            "HBD": hbd,
            "HBA": hba,
            "RotBonds": rot_bonds,
            "TPSA": round(tpsa, 2)
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccncc1"  # 吡啶
    print(f"Rule of Three 检查: {level_function(smiles)}")
