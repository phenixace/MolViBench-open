from rdkit import Chem
from rdkit.Chem import Descriptors

def level_function(mol):
    """
    计算分子的 QED 值（药物相似性评分）。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return Descriptors.qed(mol)
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CC[C@H](F)C(=O)O"
    print(f"QED: {level_function(smiles)}")