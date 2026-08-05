from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None


        pattern = Chem.MolFromSmarts("[OH]c1ccccc1")
        has_phenol = mol_obj.HasSubstructMatch(pattern)

        if not has_phenol:
            return None


        rxn = AllChem.ReactionFromSmarts('[cH:1]>>[c:1]S(=O)(=O)O')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)


        formula = rdMolDescriptors.CalcMolFormula(product)

        return {
            "has_phenol": has_phenol,
            "product": product_smiles,
            "formula": formula
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'Oc1ccccc1'
    print(f'Output: {level_function(smiles)}')
