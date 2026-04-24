from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def level_function(mol):
    """
    计算分子的TPSA（拓扑极性表面积）。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return rdMolDescriptors.CalcTPSA(mol)
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CCO"
    print(f"TPSA值: {level_function(smiles)}")