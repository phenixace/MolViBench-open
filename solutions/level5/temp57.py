import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, BRICS
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker


def level_function(active_smiles_list, library_smiles_list=None, n_generate=200, n_select=50, seed=42):
    """设计一个聚焦库（focused library）：给定靶标已知活性物，生成相似但多样的候选并评估覆盖度。"""
    try:
        np.random.seed(seed)

        active_mols = [Chem.MolFromSmiles(s) for s in active_smiles_list]
        active_mols = [m for m in active_mols if m is not None]
        if not active_mols:
            return None

        # Compute fingerprints for actives
        active_fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in active_mols]

        # Generate candidates via BRICS decomposition + recombination
        all_frags = set()
        for m in active_mols:
            frags = BRICS.BRICSDecompose(m)
            all_frags.update(frags)

        candidates = set()
        if len(all_frags) >= 2:
            try:
                builder = BRICS.BRICSBuild(list(all_frags)[:50])
                for i, prod in enumerate(builder):
                    if i >= n_generate:
                        break
                    try:
                        Chem.SanitizeMol(prod)
                        smi = Chem.MolToSmiles(prod)
                        if 100 < Descriptors.MolWt(Chem.MolFromSmiles(smi)) < 600:
                            candidates.add(smi)
                    except Exception:
                        continue
            except Exception:
                pass

        # Also add library molecules if provided
        if library_smiles_list:
            for s in library_smiles_list:
                m = Chem.MolFromSmiles(s)
                if m is not None:
                    candidates.add(Chem.MolToSmiles(m))

        if not candidates:
            return {"error": "No candidates generated"}

        cand_list = list(candidates)
        cand_mols = [Chem.MolFromSmiles(s) for s in cand_list]
        cand_fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in cand_mols if m is not None]
        valid_cands = [(s, m, fp) for s, m, fp in zip(cand_list, cand_mols, cand_fps) if m is not None]

        # Compute similarity to actives (max Tanimoto)
        scored = []
        for smi, mol, fp in valid_cands:
            max_sim = max(DataStructs.TanimotoSimilarity(fp, afp) for afp in active_fps)
            qed = Descriptors.qed(mol)
            scored.append({
                "smiles": smi,
                "max_similarity_to_actives": round(max_sim, 4),
                "qed": round(qed, 4),
                "mw": round(Descriptors.MolWt(mol), 2),
            })

        # Filter: similarity >= 0.2
        scored = [s for s in scored if s["max_similarity_to_actives"] >= 0.2]
        scored.sort(key=lambda x: x["max_similarity_to_actives"], reverse=True)

        # Diversity selection via MaxMin from top candidates
        top_pool = scored[:min(len(scored), n_select * 5)]
        if len(top_pool) <= n_select:
            selected = top_pool
        else:
            pool_mols = [Chem.MolFromSmiles(s["smiles"]) for s in top_pool]
            pool_fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in pool_mols]

            def dist_fn(i, j):
                return 1.0 - DataStructs.TanimotoSimilarity(pool_fps[i], pool_fps[j])

            picker = MaxMinPicker()
            picks = picker.LazyPick(dist_fn, len(pool_fps), n_select, seed=seed)
            selected = [top_pool[i] for i in picks]

        # Evaluate coverage: compute internal diversity
        if len(selected) >= 2:
            sel_fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s["smiles"]), 2, nBits=2048) for s in selected]
            dists = []
            for i in range(len(sel_fps)):
                for j in range(i + 1, len(sel_fps)):
                    dists.append(1.0 - DataStructs.TanimotoSimilarity(sel_fps[i], sel_fps[j]))
            internal_diversity = float(np.mean(dists))
        else:
            internal_diversity = 0.0

        avg_sim = float(np.mean([s["max_similarity_to_actives"] for s in selected])) if selected else 0.0

        return {
            "n_actives": len(active_mols),
            "n_candidates_generated": len(candidates),
            "n_selected": len(selected),
            "avg_similarity_to_actives": round(avg_sim, 4),
            "internal_diversity": round(internal_diversity, 4),
            "selected_molecules": selected,
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    actives = ["c1ccc(NC(=O)c2ccccc2)cc1", "c1ccc(NC(=O)c2ccncc2)cc1",
               "c1ccc(NC(=O)c2ccc(F)cc2)cc1"]
    result = level_function(actives, n_select=10)
    if result:
        print(f"Generated: {result['n_candidates_generated']}, Selected: {result['n_selected']}")
        print(f"Avg similarity: {result['avg_similarity_to_actives']}, Diversity: {result['internal_diversity']}")
