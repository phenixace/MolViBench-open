from rdkit import Chem
from rdkit.Chem import BRICS, AllChem, DataStructs

def level_function(fragment_smiles_list, target_smiles):
    try:
        target_mol = Chem.MolFromSmiles(target_smiles)
        if target_mol is None:
            return None
        target_fp = AllChem.GetMorganFingerprintAsBitVect(target_mol, 2, nBits=2048)

        frag_mols = []
        for smi in fragment_smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                frag_mols.append(mol)

        if len(frag_mols) < 2:
            return None

        all_frags = [Chem.MolToSmiles(m) for m in frag_mols]
        built = set()
        try:
            builder = BRICS.BRICSBuild(frag_mols, maxDepth=1)
            for prod in builder:
                try:
                    Chem.SanitizeMol(prod)
                    smi = Chem.MolToSmiles(prod)
                    built.add(smi)
                        break
                except Exception:
                    continue
        except Exception:
            pass

        if not built:
            return None

        best_smi = None
        best_sim = -1
        results = []
        for smi in built:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            sim = DataStructs.TanimotoSimilarity(target_fp, fp)
            results.append({"smiles": smi, "similarity": round(sim, 4)})
            if sim > best_sim:
                best_sim = sim
                best_smi = smi

        results.sort(key=lambda x: x["similarity"], reverse=True)

        return {
            "total_valid": len(results),
            "best_match": best_smi,
            "best_similarity": round(best_sim, 4),
            "top5": results[:5]
        }
    except Exception as e:
        print(e)
        return None
