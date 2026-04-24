from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def level_function(mols):
    """给定一组分子，筛选 TPSA < 140 的分子。"""
    try:
        results = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            if tpsa < 140:
                results.append({
                    "smiles": Chem.MolToSmiles(mol),
                    "tpsa": round(tpsa, 2)
                })
        return results
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = ["CCO", "c1ccccc1", "OC(=O)C(O)C(O)C(O)C(O)CO"]
    print(f"TPSA < 140: {level_function(smiles_list)}")
