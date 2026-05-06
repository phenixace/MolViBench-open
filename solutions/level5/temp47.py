from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

def level_function(fragment, linkers=None):
    try:
        frag_mol = Chem.MolFromSmiles(fragment)
        if frag_mol is None:
            return None

        if linkers is None:
            linkers = ["CC", "CCC", "CCO", "CCNC", "C(=O)"]

        candidates = set()

        growth_rxns = [
            '[cH:1]>>[c:1]CC',
            '[cH:1]>>[c:1]CCC',
            '[cH:1]>>[c:1]CCO',
            '[cH:1]>>[c:1]CCNC',
            '[cH:1]>>[c:1]CC(=O)O',
            '[cH:1]>>[c:1]Cc1ccccc1',
        ]

        for rxn_smarts in growth_rxns:
            rxn = AllChem.ReactionFromSmarts(rxn_smarts)
            products = rxn.RunReactants((frag_mol,))
            for product_set in products:
                for product in product_set:
                    try:
                        Chem.SanitizeMol(product)
                        smi = Chem.MolToSmiles(product)
                        candidates.add(smi)
                    except Exception:
                        continue

        filtered = []
        for smi in candidates:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            logp = Descriptors.MolLogP(mol)
            if 2 <= logp <= 4:
                filtered.append({'smiles': smi, 'logp': round(logp, 2), 'mol': mol})

        results = []
        for item in filtered:
            qed = Descriptors.qed(item['mol'])
            results.append({
                'smiles': item['smiles'],
                'logp': item['logp'],
                'qed': round(qed, 4)
            })

        results.sort(key=lambda x: x['qed'], reverse=True)
        return results[:10]
    except Exception as e:
        print(e)
        return None
