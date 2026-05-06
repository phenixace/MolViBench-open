from rdkit import Chem
from rdkit.Chem import Descriptors

def level_function(mols):
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
