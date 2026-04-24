from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


def level_function(mol):
    """提取分子的 Murcko 骨架（scaffold）。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        scaffold = MurckoScaffold.GetScaffoldForMol(mol_obj)
        return Chem.MolToSmiles(scaffold)
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccc(CC(=O)O)cc1"
    print(f"Murcko scaffold: {level_function(smiles)}")
