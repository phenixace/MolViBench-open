from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


def level_function(fragments):

    try:

        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog(params)


        growth_rxns = [
            '[cH:1]>>[c:1]C',
            '[cH:1]>>[c:1]CC',
            '[cH:1]>>[c:1]O',
            '[cH:1]>>[c:1]OC',
            '[cH:1]>>[c:1]F',
            '[cH:1]>>[c:1]N',
        ]

        candidates = set()
        for frag_smi in fragments:
            frag = Chem.MolFromSmiles(frag_smi)
            if frag is None:
                continue
            for rxn_smarts in growth_rxns:
                rxn = AllChem.ReactionFromSmarts(rxn_smarts)
                products = rxn.RunReactants((frag,))
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
            if catalog.GetFirstMatch(mol) is None:
                filtered.append(smi)


        results = []
        for smi in filtered:
            mol = Chem.MolFromSmiles(smi)
            logp = Descriptors.MolLogP(mol)
            results.append({
                'smiles': smi,
                'logp': round(logp, 2),
                'logp_distance': round(abs(logp - 2.5), 2)
            })

        results.sort(key=lambda x: x['logp_distance'])
        return results[:10]
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    frags = ['c1ccncc1', 'c1ccccc1', 'c1ccc2[nH]ccc2c1']
    result = level_function(frags)
    if result:
        for r in result[:5]:
            print(f"Output: {r['smiles']}{r['logp']}")
