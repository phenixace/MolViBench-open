from rdkit import Chem
from rdkit.Chem import GraphDescriptors


def level_function(mol):
    """计算分子的 BertzCT 复杂度指数。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        bertz = GraphDescriptors.BertzCT(mol_obj)
        return round(bertz, 4)
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccccc1"
    print(f"BertzCT: {level_function(smiles)}")
