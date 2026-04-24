from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol):
    """
    计算分子可形成的立体中心数。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return AllChem.CalcNumAtomStereoCenters(mol)
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CC[C@H](F)C(=O)O"
    print(f"可形成的立体中心数: {level_function(smiles)}")