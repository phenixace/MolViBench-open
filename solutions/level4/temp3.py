from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        pattern = Chem.MolFromSmarts('[NX3;H2]')
        has_amino = mol_obj.HasSubstructMatch(pattern)

        if not has_amino:
            return None

        rxn = AllChem.ReactionFromSmarts('[N:1]([H])[H]>>[N:1]C(=O)C')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        tpsa = rdMolDescriptors.CalcTPSA(product)

        return {
            "has_amino": has_amino,
            "product": product_smiles,
            "tpsa": round(tpsa, 4)
        }
    except Exception as e:
        print(e)
        return None
