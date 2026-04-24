from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol):
    """计算 3D 构象的能量。"""
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        if result != 0:
            return None
        ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
        if ff is None:
            ff = AllChem.UFFGetMoleculeForceField(mol)
        if ff is None:
            return None
        energy = ff.CalcEnergy()
        return energy
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "CCO"
    result = level_function(smiles)
    print(f"3D 构象能量: {result}")
