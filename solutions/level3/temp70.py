from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate

def level_function(mol1, mol2):
    try:
        m1 = Chem.MolFromSmiles(mol1)
        m2 = Chem.MolFromSmiles(mol2)
        if m1 is None or m2 is None:
            return None

        AllChem.Compute2DCoords(m1)
        AllChem.Compute2DCoords(m2)

        factory = Gobbi_Pharm2D.factory
        fp1 = Generate.Gen2DFingerprint(m1, factory)
        fp2 = Generate.Gen2DFingerprint(m2, factory)

        tanimoto = DataStructs.TanimotoSimilarity(fp1, fp2)
        dice = DataStructs.DiceSimilarity(fp1, fp2)

        return {
            "tanimoto": round(tanimoto, 4),
            "dice": round(dice, 4),
            "fp1_on_bits": len(fp1.GetOnBits()),
            "fp2_on_bits": len(fp2.GetOnBits())
        }
    except Exception as e:
        print(e)
        return None
