from rdkit import Chem
from rdkit.Chem import BRICS, AllChem, DataStructs


def level_function(fragment_smiles_list, target_smiles):
    """给定一组候选片段 → 尝试所有两两 BRICS 拼接 → 过滤掉化学上无效的分子 → 对有效分子计算与目标的相似度 → 返回最相似的。"""
    try:
        target_mol = Chem.MolFromSmiles(target_smiles)
        if target_mol is None:
            return None
        target_fp = AllChem.GetMorganFingerprintAsBitVect(target_mol, 2, nBits=2048)

        # Parse fragments
        frag_mols = []
        for smi in fragment_smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                frag_mols.append(mol)

        if len(frag_mols) < 2:
            return None

        # Try all pairwise BRICS builds
        all_frags = [Chem.MolToSmiles(m) for m in frag_mols]
        built = set()
        try:
            # Use BRICS.BRICSBuild for fragment assembly
            builder = BRICS.BRICSBuild(frag_mols, maxDepth=1)
            for prod in builder:
                try:
                    Chem.SanitizeMol(prod)
                    smi = Chem.MolToSmiles(prod)
                    built.add(smi)
                    if len(built) >= 500:  # Limit
                        break
                except Exception:
                    continue
        except Exception:
            pass

        if not built:
            return None

        # Calculate similarity to target
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


if __name__ == "__main__":
    frags = ["[1*]c1ccccc1", "[1*]CC(=O)O", "[1*]N", "[1*]CCO", "[1*]c1ccncc1"]
    target = "c1ccc(CC(=O)O)cc1"
    result = level_function(frags, target)
    if result:
        print(f"Best: {result['best_match']}, Sim: {result['best_similarity']}")
