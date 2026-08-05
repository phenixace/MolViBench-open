from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(product, reverse_reactions):

    try:
        mol = Chem.MolFromSmiles(product)
        if mol is None:
            return None
        all_reactants = {}
        for name, smarts in reverse_reactions.items():
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((mol,))
            if products:
                reactant_pairs = []
                for product_set in products:
                    pair = []
                    for reactant in product_set:
                        Chem.SanitizeMol(reactant)
                        pair.append(Chem.MolToSmiles(reactant))
                    if pair not in reactant_pairs:
                        reactant_pairs.append(pair)
                if reactant_pairs:
                    all_reactants[name] = reactant_pairs
        return all_reactants if all_reactants else None
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    product = 'CC(=O)OCC'
    reverse_reactions = {'ester_hydrolysis': '[C:1](=O)[O:2][C:3]>>[C:1](=O)[OH].[OH][C:3]', 'amide_hydrolysis': '[C:1](=O)[N:2]>>[C:1](=O)[OH].[N:2]'}
    print(f'Output: {level_function(product, reverse_reactions)}')
