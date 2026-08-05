from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        orig_logp = Descriptors.MolLogP(mol_obj)
        orig_smi = Chem.MolToSmiles(mol_obj)


        if orig_logp < 2:

            rxns = [
                ('[cH:1]>>[c:1]C', 'Add methyl group'),
                ('[cH:1]>>[c:1]F', 'Add fluorine'),
                ('[cH:1]>>[c:1]Cl', 'Add chlorine'),
            ]
        elif orig_logp > 3:

            rxns = [
                ('[cH:1]>>[c:1]O', 'Add hydroxyl group'),
                ('[cH:1]>>[c:1]N', 'Add amino group'),
                ('[CH3:1]>>[CH2:1]O', 'Alkyl oxidation'),
            ]
        else:
            return {"smiles": orig_smi, "logp": round(orig_logp, 2),
                    "message": "LogP is already within the 2-3 range"}

        candidates = []
        for smarts, desc in rxns:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((mol_obj,))
            for product_set in products:
                for product in product_set:
                    try:
                        Chem.SanitizeMol(product)
                        smi = Chem.MolToSmiles(product)
                        if smi == orig_smi:
                            continue
                        logp = Descriptors.MolLogP(product)
                        if 2 <= logp <= 3:
                            candidates.append({
                                'smiles': smi,
                                'logp': round(logp, 2),
                                'modification': desc
                            })
                    except Exception:
                        continue


        candidates.sort(key=lambda x: abs(x['logp'] - 2.5))
        return candidates[:5] if candidates else None
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccccc1'
    print(f'Output: {level_function(smiles)}')
