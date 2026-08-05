from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


def level_function(mol, diene="C=CC=C"):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None


        pattern = Chem.MolFromSmarts('[C:1]=[C:2]')
        has_double_bond = mol_obj.HasSubstructMatch(pattern)

        if not has_double_bond:
            return None


        diene_mol = Chem.MolFromSmiles(diene)
        if diene_mol is None:
            return None

        rxn = AllChem.ReactionFromSmarts(
            '[C:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2]=[C:3][C:4][C:6][C:5]1'
        )
        products = rxn.RunReactants((diene_mol, mol_obj))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)


        mol_wt = rdMolDescriptors.CalcExactMolWt(product)

        return {
            "has_double_bond": has_double_bond,
            "product": product_smiles,
            "molecular_weight": round(mol_wt, 4)
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'C=CC(=O)C'
    print(f'Output: {level_function(smiles)}')
