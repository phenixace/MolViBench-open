from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol_smi):



    try:
        mol_obj = Chem.MolFromSmiles(mol_smi)
        if mol_obj is None:
            return None



        pattern = Chem.MolFromSmarts('[NH2,NH1]')
        has_amine = mol_obj.HasSubstructMatch(pattern)

        if not has_amine:
            return None




        rxn = AllChem.ReactionFromSmarts('[*:1][NH2:2]>>[*:1][N+:2](=O)[O-]')




        products = rxn.RunReactants((mol_obj,))
        if not products:
            return {"has_amine": True, "product": "Reaction failed to run", "molecular_weight": None}


        product = products[0][0]


        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)


        mol_wt = rdMolDescriptors.CalcExactMolWt(product)

        return {
            "has_amine": True,
            "product": product_smiles,
            "molecular_weight": round(mol_wt, 4)
        }
    except Exception as e:
        print(f"Error logic: {e}")
        return None

if __name__ == '__main__':
    smiles = 'c1ccccc1CN'
    result = level_function(smiles)
    if result:
        print(f"Output: {result['product']}")
        print(f"Output: {result['molecular_weight']}")
