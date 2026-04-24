from rdkit import Chem
from rdkit.Chem import Descriptors


def level_function(mol):
    """计算分子的重原子数（非氢原子数）。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        return Descriptors.HeavyAtomCount(mol_obj)
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccccc1"
    print(f"重原子数: {level_function(smiles)}")
