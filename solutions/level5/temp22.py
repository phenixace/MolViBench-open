from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        orig_logp = Descriptors.MolLogP(mol_obj)
        orig_smi = Chem.MolToSmiles(mol_obj)


        hydro_rxns = [
            ('[cH:1]>>[c:1]C', 'Add methyl group'),
            ('[cH:1]>>[c:1]F', 'Add fluorine'),
            ('[cH:1]>>[c:1]Cl', 'Add chlorine'),
            ('[OH:1]>>[H:1]', 'Remove hydroxyl group'),
        ]

        derivatives = []
        for smarts, desc in hydro_rxns:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((mol_obj,))
            for product_set in products:
                for product in product_set:
                    try:
                        Chem.SanitizeMol(product)
                        smi = Chem.MolToSmiles(product)
                        if smi == orig_smi or smi in [d['smiles'] for d in derivatives]:
                            continue
                        new_logp = Descriptors.MolLogP(product)
                        if new_logp > orig_logp:
                            derivatives.append({
                                'smiles': smi,
                                'logp': round(new_logp, 2),
                                'modification': desc
                            })
                    except Exception:
                        continue

        derivatives.sort(key=lambda x: x['logp'], reverse=True)
        return derivatives[:10]
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccc(O)cc1'
    result = level_function(smiles)
    print(f'Output: {result}')
