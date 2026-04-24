from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
import sys
import os


def level_function(query_smiles, library_smiles, smarts_pattern=None, top_k=5):
    """给定 query 分子和分子库 → 子结构搜索找到匹配分子 → 计算命中分子的描述符（MW, LogP, QED, SA Score）→ 按 SA Score 升序排列 → 输出 Top-5。"""
    try:
        query = Chem.MolFromSmiles(query_smiles)
        if query is None:
            return None

        # Use query as substructure pattern, or provided SMARTS
        if smarts_pattern:
            pattern = Chem.MolFromSmarts(smarts_pattern)
        else:
            pattern = query

        if pattern is None:
            return None

        # SA Score function
        try:
            from rdkit.Chem import RDConfig
            sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
            import sascorer
            sa_func = sascorer.calculateScore
        except Exception:
            sa_func = lambda m: Descriptors.BertzCT(m) / 100.0

        # Substructure search
        hits = []
        for smi in library_smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            if mol.HasSubstructMatch(pattern):
                canonical = Chem.MolToSmiles(mol)
                mw = Descriptors.MolWt(mol)
                logp = Crippen.MolLogP(mol)
                qed = Descriptors.qed(mol)
                sa = sa_func(mol)
                hits.append({
                    "smiles": canonical,
                    "MW": round(mw, 2),
                    "LogP": round(logp, 4),
                    "QED": round(qed, 4),
                    "SA_Score": round(sa, 4)
                })

        # Sort by SA Score ascending (lower = easier to synthesize)
        hits.sort(key=lambda x: x["SA_Score"])

        return {
            "query": Chem.MolToSmiles(query),
            "total_library": len(library_smiles),
            "hits": len(hits),
            "top_k": hits[:top_k]
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    query = "c1ccccc1"
    library = ["c1ccc(O)cc1", "c1ccc(N)cc1", "c1ccc(NC(=O)C)cc1",
               "CCCCC", "CCO", "c1ccc(F)cc1", "c1ccc(-c2ccccc2)cc1"]
    result = level_function(query, library, top_k=3)
    if result:
        print(f"Hits: {result['hits']}")
        for h in result['top_k']:
            print(f"  {h['smiles']}: SA={h['SA_Score']}")
