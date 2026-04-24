from rdkit import Chem
from rdkit.Chem import Descriptors


def level_function(mol):
    """计算分子的 Labute 近似表面积（ASA）。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        asa = Descriptors.LabuteASA(mol_obj)
        return round(asa, 4)
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccccc1"
    print(f"Labute ASA: {level_function(smiles)}")
