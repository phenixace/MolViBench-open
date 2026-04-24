from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def level_function(mol):
    """
    计算给定分子的分子式。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return rdMolDescriptors.CalcMolFormula(mol)
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CCO"
    print(f"分子式: {level_function(smiles)}")