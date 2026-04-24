from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def level_function(mol):
    """
    计算分子旋转键数量。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return rdMolDescriptors.CalcNumRotatableBonds(mol)
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CCOc1ccccc1[O-]"
    print(f"旋转键数量: {level_function(smiles)}")