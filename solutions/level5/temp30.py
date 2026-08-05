from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol_obj)
        if num_aromatic_rings == 0:
            return None

        orig_logp = Descriptors.MolLogP(mol_obj)
        orig_smi = Chem.MolToSmiles(mol_obj)

        substituents = [
            ('[cH:1]>>[c:1]C', 'Methyl'),
            ('[cH:1]>>[c:1]O', 'Hydroxyl'),
            ('[cH:1]>>[c:1]N', 'Amino'),
            ('[cH:1]>>[c:1]F', 'Fluoro'),
            ('[cH:1]>>[c:1]Cl', 'Chloro'),
            ('[cH:1]>>[c:1]OC', 'Methoxy'),
            ('[cH:1]>>[c:1]C(F)(F)F', 'Trifluoromethyl'),
            ('[cH:1]>>[c:1][N+](=O)[O-]', 'Nitro'),
        ]

        results = []
        for smarts, name in substituents:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((mol_obj,))
            if products:
                product = products[0][0]
                try:
                    Chem.SanitizeMol(product)
                    smi = Chem.MolToSmiles(product)
                    logp = Descriptors.MolLogP(product)
                    results.append({
                        'substituent': name,
                        'smiles': smi,
                        'logp': round(logp, 2),
                        'delta_logp': round(logp - orig_logp, 2)
                    })
                except Exception:
                    continue

        results.sort(key=lambda x: x['logp'])
        return {"original_logp": round(orig_logp, 2), "derivatives": results}
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccccc1'
    result = level_function(smiles)
    if result:
        print(f"Output: {result['original_logp']}")
        for d in result['derivatives']:
            print(f"Output: {d['substituent']}{d['logp']}{d['delta_logp']}")
