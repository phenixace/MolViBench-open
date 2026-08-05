from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None


        pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H]')
        has_carboxyl = mol_obj.HasSubstructMatch(pattern)

        if not has_carboxyl:
            return None


        rxn = AllChem.ReactionFromSmarts('[C:1](=O)[OH]>>[C:1](=O)Cl')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)


        mol_wt = rdMolDescriptors.CalcExactMolWt(product)

        return {
            "has_carboxyl": has_carboxyl,
            "product": product_smiles,
            "molecular_weight": round(mol_wt, 4)
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'CC(=O)O'
    print(f'Output: {level_function(smiles)}')
