from rdkit import Chem
from rdkit.Chem import BRICS, Descriptors

def level_function(smiles_list, max_products=200):
    try:
        mols = []
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                mols.append(mol)

        if not mols:
            return None

        all_frags = set()
        for mol in mols:
            frags = BRICS.BRICSDecompose(mol)
            all_frags.update(frags)

        frag_mols = []
        for frag_smi in all_frags:
            frag_mol = Chem.MolFromSmiles(frag_smi)
            if frag_mol:
                frag_mols.append(frag_mol)

        if len(frag_mols) < 2:
            return None

        builder = BRICS.BRICSBuild(frag_mols, maxDepth=1)
        products = set()
        original_set = set(Chem.MolToSmiles(m) for m in mols)

        count = 0
        for prod in builder:
            try:
                Chem.SanitizeMol(prod)
                smi = Chem.MolToSmiles(prod)
                if smi not in original_set:
                    products.add(smi)
                count += 1
                if count >= max_products:
                    break
            except Exception:
                continue

        scored = []
        for smi in products:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                qed = Descriptors.qed(mol)
                mw = Descriptors.MolWt(mol)
                scored.append({
                    "smiles": smi,
                    "QED": round(qed, 4),
                    "MW": round(mw, 2)
                })

        scored.sort(key=lambda x: x["QED"], reverse=True)

        return {
            "input_molecules": len(mols),
            "fragments": len(all_frags),
            "new_molecules": len(products),
            "top10": scored[:10],
            "avg_QED": round(sum(s["QED"] for s in scored) / len(scored), 4) if scored else 0
        }
    except Exception as e:
        print(e)
        return None
