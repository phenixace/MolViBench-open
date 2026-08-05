from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None



        rxn_primary = AllChem.ReactionFromSmarts(
            '[NH2:1]>>[NH:1]C(=O)OC(C)(C)C'
        )

        products = rxn_primary.RunReactants((mol_obj,))
        if products:
            prod = products[0][0]
            try:
                Chem.SanitizeMol(prod)
                return Chem.MolToSmiles(prod)
            except Exception:
                pass


        rxn_secondary = AllChem.ReactionFromSmarts(
            '[NH1:1]>>[N:1]C(=O)OC(C)(C)C'
        )
        products = rxn_secondary.RunReactants((mol_obj,))
        if products:
            prod = products[0][0]
            try:
                Chem.SanitizeMol(prod)
                return Chem.MolToSmiles(prod)
            except Exception:
                pass

        return None
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccc(N)cc1'
    print(f'Output: {level_function(smiles)}')
