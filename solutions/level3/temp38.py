from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


def level_function(mol1, mol2):
    """判断两个分子是否为互变异构体。"""
    try:
        m1 = Chem.MolFromSmiles(mol1)
        m2 = Chem.MolFromSmiles(mol2)
        if m1 is None or m2 is None:
            return False

        enumerator = rdMolStandardize.TautomerEnumerator()

        canon1 = enumerator.Canonicalize(m1)
        canon2 = enumerator.Canonicalize(m2)

        canon_smi1 = Chem.MolToSmiles(canon1)
        canon_smi2 = Chem.MolToSmiles(canon2)

        if canon_smi1 == canon_smi2:
            smi1 = Chem.MolToSmiles(m1)
            smi2 = Chem.MolToSmiles(m2)
            if smi1 != smi2:
                return True
            else:
                return False

        return False
    except Exception as e:
        print(e)
        return False


if __name__ == "__main__":
    mol1 = "CC(=O)CC"
    mol2 = "CC(O)=CC"
    print(f"是否为互变异构体: {level_function(mol1, mol2)}")
