from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


def _synthetic_complexity_score(mol):

    num_rings = mol.GetRingInfo().NumRings()
    num_chiral = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    num_hetero = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() not in [1, 6])
    mw = Descriptors.MolWt(mol)
    return round(num_rings * 1.5 + num_chiral * 2.0 + num_hetero * 0.5 + mw / 200.0, 2)


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None


        halogen_rxns = [
            ('[cH:1]>>[c:1]F', 'F'),
            ('[cH:1]>>[c:1]Cl', 'Cl'),
            ('[cH:1]>>[c:1]Br', 'Br'),
        ]

        derivatives = set()
        for rxn_smarts, _ in halogen_rxns:
            rxn = AllChem.ReactionFromSmarts(rxn_smarts)
            products = rxn.RunReactants((mol_obj,))
            for product_set in products:
                for product in product_set:
                    try:
                        Chem.SanitizeMol(product)
                        smi = Chem.MolToSmiles(product)
                        derivatives.add(smi)
                    except Exception:
                        continue


        filtered = []
        for smi in derivatives:
            m = Chem.MolFromSmiles(smi)
            if m is None:
                continue
            mw = Descriptors.MolWt(m)
            if mw < 450:
                filtered.append({'smiles': smi, 'mw': round(mw, 2), 'mol': m})


        results = []
        for item in filtered:
            score = _synthetic_complexity_score(item['mol'])
            results.append({
                'smiles': item['smiles'],
                'mw': item['mw'],
                'complexity_score': score
            })

        results.sort(key=lambda x: x['complexity_score'])
        return results[:10]
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccc(CC(=O)O)cc1'
    result = level_function(smiles)
    if result:
        for r in result[:5]:
            print(f"Output: {r['smiles']}{r['mw']}{r['complexity_score']}")
