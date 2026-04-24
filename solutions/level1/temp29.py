from rdkit import Chem
from rdkit.Chem import Descriptors

def level_function(mol):
    """
    获取分子的原子价电子数分布。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return Descriptors.NumValenceElectrons(mol)
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CC[C@H](F)C(=O)O"
    print(f"原子价电子数分布: {level_function(smiles)}")