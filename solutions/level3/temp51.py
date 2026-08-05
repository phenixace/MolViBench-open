from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(product_smiles):

    try:
        product = Chem.MolFromSmiles(product_smiles)
        if product is None:
            return None


        retro_rxns = {
            "amide_bond": "[C:1](=[O:2])[NH:3]>>[C:1](=[O:2])O.[NH2:3]",
            "ester_bond": "[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])O.[OH:3][C:4]",
            "ether_bond": "[C:1][O:2][C:3]>>[C:1]O.[OH:2][C:3]",
            "C-N_bond": "[C:1][NH:2]>>[C:1]Br.[NH2:2]",
            "Suzuki": "[c:1][c:2]>>[c:1]Br.[c:2]B(O)O",
        }

        results = []
        for rxn_name, rxn_smarts in retro_rxns.items():
            try:
                rxn = AllChem.ReactionFromSmarts(rxn_smarts)
                if rxn is None:
                    continue
                prod_sets = rxn.RunReactants((product,))
                for prod_set in prod_sets:
                    reactants = []
                    valid = True
                    for mol in prod_set:
                        try:
                            Chem.SanitizeMol(mol)
                            reactants.append(Chem.MolToSmiles(mol))
                        except Exception:
                            valid = False
                            break
                    if valid and reactants:
                        results.append({
                            "reaction_type": rxn_name,
                            "reactants": reactants
                        })
            except Exception:
                continue

        return results if results else None
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccc(NC(=O)c2ccccc2)cc1'
    result = level_function(smiles)
    if result:
        for r in result:
            print(f"Output: {r['reaction_type']}{r['reactants']}")
