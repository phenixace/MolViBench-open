from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(filename):
    """从 PDB 文件中读取分子。"""
    try:
        mol = Chem.MolFromPDBFile(filename, removeHs=True)
        if mol is None:
            return None
        smiles = Chem.MolToSmiles(mol)
        return smiles
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    # First create a test PDB file
    mol = Chem.MolFromSmiles("CCO")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    Chem.MolToPDBFile(mol, "test_input.pdb")
    result = level_function("test_input.pdb")
    print(f"从 PDB 读取的分子: {result}")
