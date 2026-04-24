from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski


def level_function(mols):
    """给定一组分子，筛选符合 Lipinski Rule of Five 的分子。"""
    try:
        results = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Lipinski.NumHDonors(mol)
            hba = Lipinski.NumHAcceptors(mol)
            if mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10:
                results.append(Chem.MolToSmiles(mol))
        return results
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = ["CCO", "c1ccccc1", "CC(=O)O",
                   "C" * 50,  # 长链烷烃
                   "c1ccc(O)cc1"]
    print(f"Lipinski 符合: {level_function(smiles_list)}")
