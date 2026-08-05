from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol1, mol2):

    try:
        reaction_smarts = '[C:1][Cl:2].[OH2:3]>>[C:1][OH]'
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        halide = Chem.MolFromSmiles(mol1)
        nucleophile = Chem.MolFromSmiles(mol2)
        if halide is None or nucleophile is None:
            return None
        products = rxn.RunReactants((halide, nucleophile))
        if not products:
            return None
        result_smiles = []
        for product_set in products:
            for product in product_set:
                Chem.SanitizeMol(product)
                smi = Chem.MolToSmiles(product)
                if smi not in result_smiles:
                    result_smiles.append(smi)
        return result_smiles
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    halide = 'CCCl'
    nucleophile = 'O'
    print(f'Output: {level_function(halide, nucleophile)}')
