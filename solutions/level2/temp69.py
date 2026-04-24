from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(reactants_list1, reactants_list2, reaction_smarts):
    """给定两组反应物列表和一个反应 SMARTS，枚举所有可能的产物并返回去重后的 SMILES 列表。"""
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


if __name__ == "__main__":
    # Amide coupling: acid + amine -> amide
    rxn = "[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]"
    acids = ["CC(=O)O", "c1ccc(C(=O)O)cc1"]
    amines = ["CN", "CCN", "c1ccc(N)cc1"]
    result = level_function(acids, amines, rxn)
    print(f"产物数: {len(result) if result else 0}")
    if result:
        for smi in result:
            print(f"  {smi}")
