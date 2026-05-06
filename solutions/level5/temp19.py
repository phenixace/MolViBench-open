from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def level_function(mols):
    try:

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
