from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        pattern = Chem.MolFromSmarts('[F,Cl,Br,I]')
        has_halogen = mol_obj.HasSubstructMatch(pattern)

        if not has_halogen:
            return None

        rxn = AllChem.ReactionFromSmarts('[C:1][Cl,Br,I,F]>>[C:1]O')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        mol_wt = rdMolDescriptors.CalcExactMolWt(product)

        return {
            "has_halogen": has_halogen,
            "product": product_smiles,
            "molecular_weight": round(mol_wt, 4)
        }
    except Exception as e:
        print(e)
        return None
