from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors

def level_function(scaffold, n_modifications=10):
    try:
        scaffold_mol = Chem.MolFromSmiles(scaffold)
        if scaffold_mol is None:
            return None

        scaffold_fp = AllChem.GetMorganFingerprintAsBitVect(scaffold_mol, 2, nBits=2048)

        rxns = [
            ('[cH:1]>>[c:1]C', 'Methyl'),
            ('[cH:1]>>[c:1]CC', 'Ethyl'),
            ('[cH:1]>>[c:1]O', 'Hydroxyl'),
            ('[cH:1]>>[c:1]OC', 'Methoxy'),
            ('[cH:1]>>[c:1]F', 'Fluorine'),
            ('[cH:1]>>[c:1]Cl', 'Chlorine'),
            ('[cH:1]>>[c:1]N', 'Amino'),
            ('[cH:1]>>[c:1]NC(=O)C', 'Amide'),
            ('[cH:1]>>[c:1]C(F)(F)F', 'Trifluoromethyl'),
            ('[cH:1]>>[c:1]C#N', 'Cyano'),
        ]

        derivatives = []
        seen = set()
        for rxn_smarts, name in rxns:
            rxn = AllChem.ReactionFromSmarts(rxn_smarts)
            products = rxn.RunReactants((scaffold_mol,))
            for product_set in products:
                for product in product_set:
                    try:
                        Chem.SanitizeMol(product)
                        smi = Chem.MolToSmiles(product)
                        if smi in seen:
                            continue
                        seen.add(smi)

                        qed = Descriptors.qed(product)
                        if qed > 0.7:
                            fp = AllChem.GetMorganFingerprintAsBitVect(product, 2, nBits=2048)
                            sim = DataStructs.TanimotoSimilarity(scaffold_fp, fp)
                            derivatives.append({
                                'smiles': smi,
                                'qed': round(qed, 4),
                                'similarity_to_scaffold': round(sim, 4),
                                'modification': name
                            })
                    except Exception:
                        continue

        derivatives.sort(key=lambda x: x['qed'], reverse=True)
        return derivatives[:n_modifications]
    except Exception as e:
        print(e)
        return None
