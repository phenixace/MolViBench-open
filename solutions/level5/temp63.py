import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_squared_error


def level_function(smiles_list, activities, fp_radius=2, fp_bits=2048, seed=42):

    try:

        data = []
        for smi, act in zip(smiles_list, activities):
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            try:
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                scaffold_smi = Chem.MolToSmiles(scaffold)
            except Exception:
                scaffold_smi = ""
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, fp_radius, nBits=fp_bits)
            data.append({
                "smiles": smi,
                "activity": float(act),
                "scaffold": scaffold_smi,
                "fp": np.array(fp, dtype=int),
            })

        if len(data) < 10:
            return None



        scaffold_groups = {}
        for i, d in enumerate(data):
            sc = d["scaffold"]
            if sc not in scaffold_groups:
                scaffold_groups[sc] = []
            scaffold_groups[sc].append(i)


        scaffolds_sorted = sorted(scaffold_groups.keys(), key=lambda x: len(scaffold_groups[x]), reverse=True)


        train_idx = []
        test_idx = []
        n_total = len(data)
        np.random.seed(seed)

        for sc in scaffolds_sorted:
            indices = scaffold_groups[sc]
            if len(train_idx) < 0.8 * n_total:
                train_idx.extend(indices)
            else:
                test_idx.extend(indices)


        if len(test_idx) == 0:

            last_sc = scaffolds_sorted[-1]
            moved = scaffold_groups[last_sc]
            train_idx = [i for i in train_idx if i not in moved]
            test_idx = moved


        X_train = np.array([data[i]["fp"] for i in train_idx])
        y_train = np.array([data[i]["activity"] for i in train_idx])
        X_test = np.array([data[i]["fp"] for i in test_idx])
        y_test = np.array([data[i]["activity"] for i in test_idx])


        rf = RandomForestRegressor(n_estimators=100, random_state=seed)
        rf.fit(X_train, y_train)


        train_pred = rf.predict(X_train)
        test_pred = rf.predict(X_test)

        train_r2 = r2_score(y_train, train_pred)
        train_rmse = float(np.sqrt(mean_squared_error(y_train, train_pred)))
        test_r2 = r2_score(y_test, test_pred)
        test_rmse = float(np.sqrt(mean_squared_error(y_test, test_pred)))


        n_scaffolds = len(scaffold_groups)
        train_scaffolds = set(data[i]["scaffold"] for i in train_idx)
        test_scaffolds = set(data[i]["scaffold"] for i in test_idx)
        scaffold_overlap = len(train_scaffolds & test_scaffolds)

        return {
            "n_molecules": len(data),
            "n_scaffolds": n_scaffolds,
            "split": {
                "n_train": len(train_idx),
                "n_test": len(test_idx),
                "train_scaffolds": len(train_scaffolds),
                "test_scaffolds": len(test_scaffolds),
                "scaffold_overlap": scaffold_overlap,
            },
            "train_metrics": {
                "R2": round(train_r2, 4),
                "RMSE": round(train_rmse, 4),
            },
            "test_metrics": {
                "R2": round(test_r2, 4),
                "RMSE": round(test_rmse, 4),
            },
            "test_predictions": [
                {
                    "smiles": data[i]["smiles"],
                    "actual": round(data[i]["activity"], 4),
                    "predicted": round(float(p), 4),
                }
                for i, p in zip(test_idx, test_pred)
            ][:20],
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = ['c1ccccc1', 'c1ccc(O)cc1', 'c1ccc(N)cc1', 'c1ccc(F)cc1', 'c1ccc(Cl)cc1', 'c1ccc(Br)cc1', 'c1ccc(C)cc1', 'c1ccc(OC)cc1', 'c1ccc(C(=O)O)cc1', 'c1ccc(C#N)cc1', 'c1ccncc1', 'c1ccc2ccccc2c1', 'c1ccc(CC)cc1', 'c1ccc(CO)cc1', 'CC(=O)O', 'CCCO', 'CCO', 'CC(C)O', 'CC(=O)N', 'CCC(=O)O']
    acts = [1.0, 1.5, 2.0, 1.2, 1.8, 2.5, 0.8, 1.4, 3.0, 2.2, 1.1, 1.6, 0.9, 1.3, 0.5, 0.6, 0.4, 0.7, 1.0, 0.8]
    result = level_function(smiles, acts)
    if result:
        print(f"Output: {result['train_metrics']['R2']}{result['train_metrics']['RMSE']}")
        print(f"Output: {result['test_metrics']['R2']}{result['test_metrics']['RMSE']}")
        print(f"Output: {result['split']}")
