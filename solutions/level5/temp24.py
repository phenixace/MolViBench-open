from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        orig_tpsa = rdMolDescriptors.CalcTPSA(mol_obj)
        orig_smi = Chem.MolToSmiles(mol_obj)

        if orig_tpsa < 100:
            return {"smiles": orig_smi, "tpsa": round(orig_tpsa, 2),
                    "message": "TPSA already < 100"}

        rxns = [
            ('[OH:1]>>[OCH3:1]', 'Hydroxyl methylation'),
            ('[NH2:1]>>[NHC:1]', 'Amino alkylation'),
            ('[C:1](=O)[OH]>>[C:1](=O)OC', 'Esterification'),
        ]

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
                        new_tpsa = rdMolDescriptors.CalcTPSA(product)
                        if new_tpsa < orig_tpsa:
                            candidates.append({
                                'smiles': smi,
                                'tpsa': round(new_tpsa, 2),
                                'modification': desc
                            })
                    except Exception:
                        continue

        candidates.sort(key=lambda x: x['tpsa'])
        return candidates[:5] if candidates else None
    except Exception as e:
        print(e)
        return None
