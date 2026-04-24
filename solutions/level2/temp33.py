from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(filename):
    """从 MOL 文件中读取分子。"""
    try:
        mol = Chem.MolFromMolFile(filename, removeHs=True)
        if mol is None:
            return None
        smiles = Chem.MolToSmiles(mol)
        return smiles
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    # First create a test MOL file
    mol = Chem.MolFromSmiles("CCO")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    Chem.MolToMolFile(mol, "test_input.mol")
    result = level_function("test_input.mol")
    print(f"从 MOL 读取的分子: {result}")
