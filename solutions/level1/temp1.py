from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def level_function(mol):
    """计算给定分子的分子量。"""
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return rdMolDescriptors.CalcExactMolWt(mol)
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CCO"
    print(f"分子量: {level_function(smiles)}")