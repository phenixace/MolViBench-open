from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        pattern = Chem.MolFromSmarts('[C:1]=[C:2]')
        has_alkene = mol_obj.HasSubstructMatch(pattern)

        if not has_alkene:
            return None

        rxn = AllChem.ReactionFromSmarts('[C:1]=[C:2]>>[C:1]1O[C:2]1')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        tpsa = rdMolDescriptors.CalcTPSA(product)

        return {
            "has_alkene": has_alkene,
            "product": product_smiles,
            "tpsa": round(tpsa, 4)
        }
    except Exception as e:
        print(e)
        return None
