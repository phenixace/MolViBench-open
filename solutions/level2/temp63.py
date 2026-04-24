from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
import pandas as pd
import io


def level_function(sdf_content):
    """从 SDF 文件批量读取分子并计算所有描述符，输出为 DataFrame。"""
    try:
        supplier = Chem.SDMolSupplier()
        supplier.SetData(sdf_content)

        rows = []
        for mol in supplier:
            if mol is None:
                continue
            smi = Chem.MolToSmiles(mol)
            desc = {
                "SMILES": smi,
                "MW": round(Descriptors.MolWt(mol), 2),
                "LogP": round(Descriptors.MolLogP(mol), 4),
                "TPSA": round(Descriptors.TPSA(mol), 2),
                "HBD": Descriptors.NumHDonors(mol),
                "HBA": Descriptors.NumHAcceptors(mol),
                "RotBonds": Descriptors.NumRotatableBonds(mol),
                "RingCount": Descriptors.RingCount(mol),
                "HeavyAtomCount": Descriptors.HeavyAtomCount(mol),
                "QED": round(Descriptors.qed(mol), 4),
            }
            rows.append(desc)

        df = pd.DataFrame(rows)
        return df
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    # Create a small test SDF in memory
    mols = [Chem.MolFromSmiles(s) for s in ["CCO", "c1ccccc1", "CC(=O)O"]]
    sdf_str = ""
    for m in mols:
        if m:
            AllChem.Compute2DCoords(m)
            sdf_str += Chem.MolToMolBlock(m) + "$$$$\n"
    result = level_function(sdf_str)
    print(result)
