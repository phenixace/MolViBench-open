from rdkit import Chem
from rdkit.Chem import Descriptors


def level_function(mols):
    """给定一组分子，筛选 QED > 0.6 的分子。"""
    try:
        results = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            qed = Descriptors.qed(mol)
            if qed > 0.6:
                results.append({
                    "smiles": Chem.MolToSmiles(mol),
                    "qed": round(qed, 4)
                })
        return results
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]
    print(f"QED > 0.6: {level_function(smiles_list)}")
