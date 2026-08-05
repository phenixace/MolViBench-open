from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol1, mol2):

    try:
        reaction_smarts = '[OH:1].[C:2][Cl,Br,I]>>[O:1][C:2]'
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        alcohol = Chem.MolFromSmiles(mol1)
        halide = Chem.MolFromSmiles(mol2)
        if alcohol is None or halide is None:
            return None
        products = rxn.RunReactants((alcohol, halide))
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
    alcohol = 'CCO'
    halide = 'CCBr'
    print(f'Output: {level_function(alcohol, halide)}')
