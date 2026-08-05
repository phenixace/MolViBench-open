from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None


        rxn = AllChem.ReactionFromSmarts(
            '[N:1]C(=O)OC(C)(C)C>>[NH2:1]'
        )

        products = rxn.RunReactants((mol_obj,))
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
    smiles = 'CC(C)(C)OC(=O)Nc1ccccc1'
    print(f'Output: {level_function(smiles)}')
