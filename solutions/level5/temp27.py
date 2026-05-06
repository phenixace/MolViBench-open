from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        orig_mw = Descriptors.MolWt(mol_obj)
        orig_smi = Chem.MolToSmiles(mol_obj)

        growth_rxns = [
            ('[cH:1]>>[c:1]C', 'Add methyl'),
            ('[cH:1]>>[c:1]OC', 'Add methoxy'),
            ('[cH:1]>>[c:1]F', 'Add fluorine'),
            ('[cH:1]>>[c:1]Cl', 'Add chlorine'),
            ('[cH:1]>>[c:1]CC', 'Add ethyl'),
            ('[NH2:1]>>[NH:1]C(=O)C', 'Acylation'),
        ]

        derivatives = []
        for smarts, desc in growth_rxns:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((mol_obj,))
            for product_set in products:
                for product in product_set:
                    try:
                        Chem.SanitizeMol(product)
                        smi = Chem.MolToSmiles(product)
                        if smi == orig_smi or smi in [d['smiles'] for d in derivatives]:
                            continue
                        new_mw = Descriptors.MolWt(product)
                        qed = Descriptors.qed(product)
                        if new_mw > orig_mw and qed > 0.5:
                            derivatives.append({
                                'smiles': smi,
                                'mw': round(new_mw, 2),
                                'qed': round(qed, 4),
                                'modification': desc
                            })
                    except Exception:
                        continue

        derivatives.sort(key=lambda x: x['qed'], reverse=True)
        return derivatives[:10]
    except Exception as e:
        print(e)
        return None
