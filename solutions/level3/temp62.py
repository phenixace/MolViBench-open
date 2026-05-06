from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(aryl_halide, amine):
    try:
        mol1 = Chem.MolFromSmiles(aryl_halide)
        mol2 = Chem.MolFromSmiles(amine)
        if mol1 is None or mol2 is None:
            return None

        rxn_primary = AllChem.ReactionFromSmarts(
            "[c:1][Cl,Br,I].[NH2:2]>>[c:1][NH:2]"
        )
        products = rxn_primary.RunReactants((mol1, mol2))

        if not products:
            rxn_secondary = AllChem.ReactionFromSmarts(
                "[c:1][Cl,Br,I].[NH1:2]>>[c:1][N:2]"
            )
            products = rxn_secondary.RunReactants((mol1, mol2))

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
