from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def level_function(mols):
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
