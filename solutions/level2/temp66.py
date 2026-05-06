from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

def level_function(query_smiles, library_smiles, top_k=5, radius=2, nBits=2048):
    try:
        query_mol = Chem.MolFromSmiles(query_smiles)
        if query_mol is None:
            return None

        query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, radius, nBits=nBits)

        results = []
        for smi in library_smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
            sim = DataStructs.TanimotoSimilarity(query_fp, fp)
            results.append({
                "smiles": Chem.MolToSmiles(mol),
                "similarity": round(sim, 4)
            })

        results.sort(key=lambda x: x["similarity"], reverse=True)
        return results[:top_k]
    except Exception as e:
        print(e)
        return None
