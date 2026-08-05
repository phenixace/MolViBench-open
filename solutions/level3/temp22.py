from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(reactants, reaction_smarts):

    try:
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        mols = []
        for smi in reactants:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                return None
            mols.append(mol)
        products = rxn.RunReactants(tuple(mols))
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
    reactants = ['CC(=O)O', 'CCN']
    reaction_smarts = '[C:1](=O)[OH].[N:2]>>[C:1](=O)[N:2]'
    print(f'Output: {level_function(reactants, reaction_smarts)}')
