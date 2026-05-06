import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score

def level_function(smiles_list, labels, new_smiles_list=None, fp_radius=2, fp_bits=2048, seed=42):
    try:
        mols = [Chem.MolFromSmiles(s) for s in smiles_list]
        if any(m is None for m in mols):
            return None

        X = np.array([
            list(AllChem.GetMorganFingerprintAsBitVect(m, fp_radius, nBits=fp_bits))
            for m in mols
        ], dtype=int)
        y = np.array(labels, dtype=int)

        clf = RandomForestClassifier(n_estimators=100, random_state=seed)

        n_cv = min(5, len(y))
        cv_acc = cross_val_score(clf, X, y, cv=n_cv, scoring='accuracy')
        cv_f1 = cross_val_score(clf, X, y, cv=n_cv, scoring='f1_macro')

        clf.fit(X, y)
        train_pred = clf.predict(X)
        train_acc = accuracy_score(y, train_pred)
        train_prec = precision_score(y, train_pred, average='macro', zero_division=0)
        train_rec = recall_score(y, train_pred, average='macro', zero_division=0)
        train_f1 = f1_score(y, train_pred, average='macro', zero_division=0)

        result = {
            "n_molecules": len(smiles_list),
            "n_active": int(np.sum(y == 1)),
            "n_inactive": int(np.sum(y == 0)),
            "cv_accuracy_mean": round(float(np.mean(cv_acc)), 4),
            "cv_accuracy_std": round(float(np.std(cv_acc)), 4),
            "cv_f1_mean": round(float(np.mean(cv_f1)), 4),
            "cv_f1_std": round(float(np.std(cv_f1)), 4),
            "train_accuracy": round(train_acc, 4),
            "train_precision": round(train_prec, 4),
            "train_recall": round(train_rec, 4),
            "train_f1": round(train_f1, 4),
        }

        if new_smiles_list:
            new_mols = [Chem.MolFromSmiles(s) for s in new_smiles_list]
            valid = [(i, m) for i, m in enumerate(new_mols) if m is not None]
            if valid:
                idxs, ms = zip(*valid)
                X_new = np.array([
                    list(AllChem.GetMorganFingerprintAsBitVect(m, fp_radius, nBits=fp_bits))
                    for m in ms
                ], dtype=int)
                preds = clf.predict(X_new)
                probas = clf.predict_proba(X_new)
                result["predictions"] = [
                    {
                        "smiles": new_smiles_list[i],
                        "predicted_label": int(p),
                        "probability_active": round(float(prob[1]) if prob.shape[0] > 1 else float(prob[0]), 4)
                    }
                    for i, p, prob in zip(idxs, preds, probas)
                ]

        return result
    except Exception as e:
        print(e)
        return None
