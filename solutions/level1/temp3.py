from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def level_function(mol):
    """
    计算分子的LogP值。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CCO"
    print(f"LogP值: {level_function(smiles)}")