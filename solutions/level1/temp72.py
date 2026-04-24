from rdkit import Chem
from rdkit.Chem import Descriptors


def level_function(mol):
    """判断分子是否符合 Veber 规则（旋转键 ≤10 且 TPSA ≤140）。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        rot_bonds = Descriptors.NumRotatableBonds(mol_obj)
        tpsa = Descriptors.TPSA(mol_obj)
        passes = rot_bonds <= 10 and tpsa <= 140
        return {
            "RotatableBonds": rot_bonds,
            "TPSA": round(tpsa, 2),
            "passes_Veber": passes
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    print(f"Veber规则: {level_function(smiles)}")
