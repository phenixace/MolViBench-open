from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None


        pattern = Chem.MolFromSmarts('[C:1]#[C:2]')
        has_alkyne = mol_obj.HasSubstructMatch(pattern)

        if not has_alkyne:
            return None


        rxn = AllChem.ReactionFromSmarts('[C:1]#[C:2]>>[C:1](Br)=[C:2]Br')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)


        formula = rdMolDescriptors.CalcMolFormula(product)

        return {
            "has_alkyne": has_alkyne,
            "product": product_smiles,
            "formula": formula
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'CC#C'
    print(f'Output: {level_function(smiles)}')
