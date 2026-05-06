import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score
from sklearn.metrics import mean_squared_error, r2_score

def level_function(smiles_list, activities, new_smiles_list=None, fp_radius=2, fp_bits=2048, seed=42):
    try:
        mols = [Chem.MolFromSmiles(s) for s in smiles_list]
        if any(m is None for m in mols):
            return None

        fps = []
        descs = []
        for m in mols:
            fp = AllChem.GetMorganFingerprintAsBitVect(m, fp_radius, nBits=fp_bits)
            fps.append(np.array(fp, dtype=int))
            descs.append([
                Descriptors.MolWt(m),
                Descriptors.MolLogP(m),
                Descriptors.TPSA(m),
                Descriptors.NumHDonors(m),
                Descriptors.NumHAcceptors(m),
                Descriptors.NumRotatableBonds(m),
                Descriptors.FractionCSP3(m),
            ])

        X_fp = np.array(fps)
        X_desc = np.array(descs)
        X = np.hstack([X_fp, X_desc])
        y = np.array(activities, dtype=float)

        rf = RandomForestRegressor(n_estimators=100, random_state=seed)

        cv_scores = cross_val_score(rf, X, y, cv=min(5, len(y)), scoring='r2')

        rf.fit(X, y)
        train_pred = rf.predict(X)
        train_r2 = r2_score(y, train_pred)
        train_rmse = float(np.sqrt(mean_squared_error(y, train_pred)))

        result = {
            "n_molecules": len(smiles_list),
            "n_features": X.shape[1],
            "cv_r2_mean": round(float(np.mean(cv_scores)), 4),
            "cv_r2_std": round(float(np.std(cv_scores)), 4),
            "train_r2": round(train_r2, 4),
            "train_rmse": round(train_rmse, 4),
        }

        if new_smiles_list:
            new_mols = [Chem.MolFromSmiles(s) for s in new_smiles_list]
            new_fps = []
            new_descs = []
            valid_idx = []
            for i, m in enumerate(new_mols):
                if m is None:
                    continue
                valid_idx.append(i)
                fp = AllChem.GetMorganFingerprintAsBitVect(m, fp_radius, nBits=fp_bits)
                new_fps.append(np.array(fp, dtype=int))
                new_descs.append([
                    Descriptors.MolWt(m),
                    Descriptors.MolLogP(m),
                    Descriptors.TPSA(m),
                    Descriptors.NumHDonors(m),
                    Descriptors.NumHAcceptors(m),
                    Descriptors.NumRotatableBonds(m),
                    Descriptors.FractionCSP3(m),
                ])

            if new_fps:
                X_new = np.hstack([np.array(new_fps), np.array(new_descs)])
                preds = rf.predict(X_new)
                result["predictions"] = [
                    {"smiles": new_smiles_list[i], "predicted_activity": round(float(p), 4)}
                    for i, p in zip(valid_idx, preds)
                ]

        return result
    except Exception as e:
        print(e)
        return None
