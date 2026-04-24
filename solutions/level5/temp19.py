from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


def level_function(mols):
    """给定一组分子，筛选符合 CNS 药物分布规则的分子。"""
    try:
        # CNS 药物规则 (simplified):
        # MW ≤ 400
        # LogP 1~5
        # HBD ≤ 3
        # HBA ≤ 7
        # TPSA ≤ 90
        # RotBonds ≤ 8

        results = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = rdMolDescriptors.CalcNumHBD(mol)
            hba = rdMolDescriptors.CalcNumHBA(mol)
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            rot = rdMolDescriptors.CalcNumRotatableBonds(mol)

            if (mw <= 400 and 1 <= logp <= 5 and hbd <= 3 and
                    hba <= 7 and tpsa <= 90 and rot <= 8):
                results.append({
                    "smiles": Chem.MolToSmiles(mol),
                    "MW": round(mw, 2),
                    "LogP": round(logp, 2),
                    "TPSA": round(tpsa, 2)
                })
        return results
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = ["c1ccccc1", "CCO", "c1ccc(F)cc1NC(=O)C"]
    print(f"CNS 符合: {level_function(smiles_list)}")
