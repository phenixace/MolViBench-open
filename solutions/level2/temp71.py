from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate


def level_function(mol):
    """计算分子的 2D 药效团指纹（Pharm2D fingerprint）。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Need 2D coordinates for Pharm2D
        AllChem.Compute2DCoords(mol_obj)

        # Generate Pharm2D fingerprint using Gobbi factory
        factory = Gobbi_Pharm2D.factory
        fp = Generate.Gen2DFingerprint(mol_obj, factory)

        on_bits = list(fp.GetOnBits())
        return {
            "num_bits": fp.GetNumBits(),
            "num_on_bits": len(on_bits),
            "on_bits": on_bits
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccc(NC(=O)c2ccccc2)cc1"
    result = level_function(smiles)
    if result:
        print(f"Pharm2D: {result['num_bits']} bits, {result['num_on_bits']} on")
