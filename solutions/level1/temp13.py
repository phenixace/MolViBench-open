from rdkit import Chem
from rdkit.Chem import rdinchi

def level_function(mol):
    """
    输出分子的 InChI 表示。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return rdinchi.MolToInchi(mol)[0]
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CCO"
    print(f"Inchi of the SMILES: {level_function(smiles)}")