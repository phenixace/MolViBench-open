from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(aryl_halide, aryl_boronic_acid):

    try:
        mol1 = Chem.MolFromSmiles(aryl_halide)
        mol2 = Chem.MolFromSmiles(aryl_boronic_acid)
        if mol1 is None or mol2 is None:
            return None


        rxn_smarts = "[c:1][Cl,Br,I].[c:2]B(O)O>>[c:1][c:2]"
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)

        products = rxn.RunReactants((mol1, mol2))
        if not products:
            return None

        results = set()
        for prod_set in products:
            for prod in prod_set:
                try:
                    Chem.SanitizeMol(prod)
                    results.add(Chem.MolToSmiles(prod))
                except Exception:
                    pass

        return sorted(list(results)) if results else None
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    aryl_x = 'c1ccc(Br)cc1'
    aryl_b = 'c1ccc(B(O)O)cc1'
    result = level_function(aryl_x, aryl_b)
    print(f'Output: {result}')
