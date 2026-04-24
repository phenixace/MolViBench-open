from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def level_function(mol):
    """判断分子是否含有芳香环。"""
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        return num_aromatic_rings > 0
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CCO"
    print(f"是否含有芳香环: {level_function(smiles)}")