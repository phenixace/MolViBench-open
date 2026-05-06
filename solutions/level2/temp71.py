from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        AllChem.Compute2DCoords(mol_obj)

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
