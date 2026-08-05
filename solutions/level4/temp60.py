from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol, reaction_smarts_list):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        results = []
        skipped = []

        for i, rxn_smarts in enumerate(reaction_smarts_list):
            try:
                rxn = AllChem.ReactionFromSmarts(rxn_smarts)
                if rxn is None:
                    skipped.append({"index": i, "smarts": rxn_smarts, "reason": "invalid_SMARTS"})
                    continue

                products = rxn.RunReactants((mol_obj,))
                if not products:
                    skipped.append({"index": i, "smarts": rxn_smarts, "reason": "no_products"})
                    continue

                prod_smiles = set()
                for prod_set in products:
                    for prod in prod_set:
                        try:
                            Chem.SanitizeMol(prod)
                            prod_smiles.add(Chem.MolToSmiles(prod))
                        except Exception:
                            pass

                if prod_smiles:
                    results.append({
                        "index": i,
                        "smarts": rxn_smarts,
                        "products": sorted(list(prod_smiles))
                    })
                else:
                    skipped.append({"index": i, "smarts": rxn_smarts, "reason": "sanitize_failed"})
            except Exception as ex:
                skipped.append({"index": i, "smarts": rxn_smarts, "reason": str(ex)})

        return {
            "input_smiles": Chem.MolToSmiles(mol_obj),
            "successful_reactions": len(results),
            "skipped_reactions": len(skipped),
            "products": results,
            "skipped": skipped
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccc(O)cc1'
    rxns = ['[OH:1]>>[Cl:1]', '[NH2:1]>>[N:1]C(=O)C', '[cH:1]>>[c:1]C']
    result = level_function(smiles, rxns)
    if result:
        print(f"Output: {result['successful_reactions']}{result['skipped_reactions']}")
