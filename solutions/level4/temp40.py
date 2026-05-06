from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol_obj)
        has_aromatic_ring = num_aromatic_rings > 0

        if not has_aromatic_ring:
            return None

        rxn = AllChem.ReactionFromSmarts('[cH:1]>>[c:1][N+](=O)[O-]')
        products = rxn.RunReactants((mol_obj,))
        if not products or len(products) < 1:
            return None

        first_product = products[0][0]
        Chem.SanitizeMol(first_product)

        products2 = rxn.RunReactants((first_product,))
        if not products2:
            return None

        final_product = products2[0][0]
        Chem.SanitizeMol(final_product)
        product_smiles = Chem.MolToSmiles(final_product)

        formula = rdMolDescriptors.CalcMolFormula(final_product)

        return {
            "has_aromatic_ring": has_aromatic_ring,
            "product": product_smiles,
            "formula": formula
        }
    except Exception as e:
        print(e)
        return None
