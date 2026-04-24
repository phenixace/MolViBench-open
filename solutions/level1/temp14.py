from rdkit import Chem
from rdkit.Chem import Draw

def level_function(mol):
    """
    将分子可视化为 2D SVG。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return Draw.MolToImage(mol)
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CCO"
    print(f"2D SVG: {level_function(smiles)}")