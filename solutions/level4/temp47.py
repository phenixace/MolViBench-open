from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        pattern = Chem.MolFromSmarts('[NX3;H2,H1]')
        has_amine = mol_obj.HasSubstructMatch(pattern)

        if not has_amine:
            return None

        rxn = AllChem.ReactionFromSmarts('[NH2:1]>>[NH:1]C(=O)C')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            rxn = AllChem.ReactionFromSmarts('[NH1:1]>>[N:1]C(=O)C')
            products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        logp = rdMolDescriptors.CalcCrippenDescriptors(product)[0]

        return {
            "has_amine": has_amine,
            "product": product_smiles,
            "logp": round(logp, 4)
        }
    except Exception as e:
        print(e)
        return None
