from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
import sys
import os


def level_function(query_smiles, library_smiles, smarts_pattern=None, top_k=5):

    try:
        query = Chem.MolFromSmiles(query_smiles)
        if query is None:
            return None


        if smarts_pattern:
            pattern = Chem.MolFromSmarts(smarts_pattern)
        else:
            pattern = query

        if pattern is None:
            return None


        try:
            from rdkit.Chem import RDConfig
            sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
            import sascorer
            sa_func = sascorer.calculateScore
        except Exception:
            sa_func = lambda m: Descriptors.BertzCT(m) / 100.0


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


if __name__ == '__main__':
    query = 'c1ccccc1'
    library = ['c1ccc(O)cc1', 'c1ccc(N)cc1', 'c1ccc(NC(=O)C)cc1', 'CCCCC', 'CCO', 'c1ccc(F)cc1', 'c1ccc(-c2ccccc2)cc1']
    result = level_function(query, library, top_k=3)
    if result:
        print(f"Output: {result['hits']}")
        for h in result['top_k']:
            print(f"Output: {h['smiles']}{h['SA_Score']}")
