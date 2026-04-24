from rdkit import Chem
from rdkit.Chem import Descriptors


def level_function(mols):
    """给定一组分子，筛选分子量 < 500 Da 的分子。"""
    try:
        results = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            mw = Descriptors.MolWt(mol)
            if mw < 500:
                results.append({
                    "smiles": Chem.MolToSmiles(mol),
                    "mw": round(mw, 2)
                })
        return results
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = ["CCO", "c1ccccc1", "CC(=O)O"]
    print(f"MW < 500: {level_function(smiles_list)}")
