from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        orig_smi = Chem.MolToSmiles(mol_obj)

        rxns = [
            '[cH:1]>>[c:1]C',
            '[cH:1]>>[c:1]O',
            '[cH:1]>>[c:1]N',
            '[cH:1]>>[c:1]F',
            '[cH:1]>>[c:1]Cl',
            '[cH:1]>>[c:1]OC',
            '[cH:1]>>[c:1]C(F)(F)F',
            '[cH:1]>>[c:1]C#N',
            '[cH:1]>>[c:1]CC',
            '[cH:1]>>[c:1]NC(=O)C',
        ]

        candidates = set()
        candidates.add(orig_smi)

        for rxn_smarts in rxns:
            rxn = AllChem.ReactionFromSmarts(rxn_smarts)
            products = rxn.RunReactants((mol_obj,))
            for product_set in products:
                for product in product_set:
                    try:
                        Chem.SanitizeMol(product)
                        smi = Chem.MolToSmiles(product)
                        candidates.add(smi)

                        products2 = rxn.RunReactants((product,))
                        for ps2 in products2:
                            for p2 in ps2:
                                try:
                                    Chem.SanitizeMol(p2)
                                    smi2 = Chem.MolToSmiles(p2)
                                    candidates.add(smi2)
                                except Exception:
                                    continue
                    except Exception:
                        continue

        scored = []
        for smi in candidates:
            m = Chem.MolFromSmiles(smi)
            if m is None:
                continue
            qed = Descriptors.qed(m)
            logp = Descriptors.MolLogP(m)
            tpsa = rdMolDescriptors.CalcTPSA(m)

            tpsa_score = 1.0 if tpsa < 120 else max(0, 1.0 - (tpsa - 120) / 60.0)
            total_score = qed * 0.4 + logp_score * 0.3 + tpsa_score * 0.3

            scored.append({
                'smiles': smi,
                'qed': round(qed, 4),
                'logp': round(logp, 2),
                'tpsa': round(tpsa, 2),
                'total_score': round(total_score, 4)
            })

        scored.sort(key=lambda x: x['total_score'], reverse=True)

        return {
            'num_candidates': len(candidates),
            'best_lead': scored[0] if scored else None,
            'top5': scored[:5]
        }
    except Exception as e:
        print(e)
        return None
