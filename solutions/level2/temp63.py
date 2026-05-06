from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
import pandas as pd
import io

def level_function(sdf_content):
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
