from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        pattern = Chem.MolFromSmarts('[CX3H1](=O)')
        has_aldehyde = mol_obj.HasSubstructMatch(pattern)

        if not has_aldehyde:
            return None

        rxn = AllChem.ReactionFromSmarts('[CH:1]=O>>[C:1](=O)O')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        tpsa = rdMolDescriptors.CalcTPSA(product)

        return {
            "has_aldehyde": has_aldehyde,
            "product": product_smiles,
            "tpsa": round(tpsa, 4)
        }
    except Exception as e:
        print(e)
        return None
