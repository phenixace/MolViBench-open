from rdkit import Chem
from rdkit.Chem import Descriptors


def level_function(mols, logp_range=(0, 5), mw_max=500):
    """给定一组分子，找出在 LogP 和 分子量 上同时满足条件的分子。"""
    try:
        results = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            logp = Descriptors.MolLogP(mol)
            mw = Descriptors.MolWt(mol)

            if logp_range[0] <= logp <= logp_range[1] and mw < mw_max:
                results.append({
                    'smiles': Chem.MolToSmiles(mol),
                    'logp': round(logp, 2),
                    'mw': round(mw, 2)
                })
        return results
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "CCCCCCCCCCCCCCCCCCCC"]
    result = level_function(smiles_list)
    print(f"LogP 和 MW 满足条件: {result}")
