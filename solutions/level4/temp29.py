from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None


        pattern = Chem.MolFromSmarts('C(=O)O[#6]')
        has_ester = mol_obj.HasSubstructMatch(pattern)

        if not has_ester:
            return None


        rxn = AllChem.ReactionFromSmarts('[C:1](=O)[O][#6]>>[C:1](=O)N')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)


        mol_formula = rdMolDescriptors.CalcMolFormula(product)

        return {
            "has_ester": has_ester,
            "product": product_smiles,
            "molecular_formula": mol_formula
        }
    except Exception as e:
        print(e)
        return None

if __name__ == '__main__':
    smiles = 'CC(=O)OCC'
    print(f'Output: {level_function(smiles)}')
