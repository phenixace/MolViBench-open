from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

def level_function(smiles_list):
    try:
        valid = {}
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                canonical = Chem.MolToSmiles(mol)
                if canonical not in valid:
                    valid[canonical] = mol

        if not valid:
            return None

        rows = []
        for smi, mol in valid.items():
            rows.append({
                "SMILES": smi,
                "MW": round(Descriptors.MolWt(mol), 2),
                "LogP": round(Descriptors.MolLogP(mol), 4),
                "TPSA": round(Descriptors.TPSA(mol), 2),
                "QED": round(Descriptors.qed(mol), 4)
            })

        df = pd.DataFrame(rows)
        df = df.sort_values("QED", ascending=False).reset_index(drop=True)

        return {
            "total_input": len(smiles_list),
            "valid_unique": len(valid),
            "dataframe": df.to_dict(orient="records"),
            "csv_string": df.to_csv(index=False)
        }
    except Exception as e:
        print(e)
        return None
