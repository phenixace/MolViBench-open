from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(reactants_list1, reactants_list2, reaction_smarts):
    try:
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        if rxn is None:
            return None

        products_set = set()
        for smi1 in reactants_list1:
            mol1 = Chem.MolFromSmiles(smi1)
            if mol1 is None:
                continue
            for smi2 in reactants_list2:
                mol2 = Chem.MolFromSmiles(smi2)
                if mol2 is None:
                    continue
                try:
                    prod_sets = rxn.RunReactants((mol1, mol2))
                    for prod_set in prod_sets:
                        for prod in prod_set:
                            try:
                                Chem.SanitizeMol(prod)
                                products_set.add(Chem.MolToSmiles(prod))
                            except Exception:
                                pass
                except Exception:
                    continue

        return sorted(list(products_set))
    except Exception as e:
        print(e)
        return None
