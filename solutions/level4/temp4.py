from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, RWMol

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        pattern = Chem.MolFromSmarts('[F,Cl,Br,I]')
        has_halogen = mol_obj.HasSubstructMatch(pattern)

        if not has_halogen:
            return None

        rxn = AllChem.ReactionFromSmarts('[C:1][F,Cl,Br,I]>>[C:1][H]')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        logp = rdMolDescriptors.CalcCrippenDescriptors(product)[0]

        return {
            "has_halogen": has_halogen,
            "product": product_smiles,
            "logp": round(logp, 4)
        }
    except Exception as e:
        print(e)
        return None
