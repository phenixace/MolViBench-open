from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol):
    """计算分子中每个原子的 Gasteiger 偏电荷。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        AllChem.ComputeGasteigerCharges(mol_obj)
        charges = {}
        for atom in mol_obj.GetAtoms():
            charge = float(atom.GetProp('_GasteigerCharge'))
            charges[atom.GetIdx()] = round(charge, 4)
        return charges
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "CCO"
    print(f"Gasteiger 偏电荷: {level_function(smiles)}")
