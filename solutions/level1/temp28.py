from rdkit import Chem
from rdkit.Chem import Descriptors

def level_function(mol):
    """
    输出分子的摩尔折射率。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return Descriptors.MolMR(mol)
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CC[C@H](F)C(=O)O"
    print(f"MR: {level_function(smiles)}")