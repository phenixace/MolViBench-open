from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol):
    """生成一个分子的 3D 可视化模型。"""
    try:
        mol = Chem.MolFromSmiles(mol)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        molblock = Chem.MolToMolBlock(mol)
        return molblock
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    result = level_function("c1ccccc1CCO")
    print(f"3D MolBlock:\n{result}")
