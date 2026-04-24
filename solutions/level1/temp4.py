from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def level_function(mol):
    """计算分子的氢键供体数。"""
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return rdMolDescriptors.CalcNumHBD(mol)
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CCO"
    print(f"氢键供体数: {level_function(smiles)}")