from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol_smi):



    try:
        mol_obj = Chem.MolFromSmiles(mol_smi)
        if mol_obj is None:
            return None



        pattern = Chem.MolFromSmarts('[OX2H]')
        has_hydroxyl = mol_obj.HasSubstructMatch(pattern)

        if not has_hydroxyl:
            return None




        rxn_smarts = '[*:1][OH:2]>>[*:1][O:2]C(=O)C'
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)

        products = rxn.RunReactants((mol_obj,))

        if not products:
            return {"has_hydroxyl": True, "product": "Esterification failed", "molecular_weight": None}


        product = products[0][0]


        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)


        mol_wt = rdMolDescriptors.CalcExactMolWt(product)

        return {
            "has_hydroxyl": True,
            "product": product_smiles,
            "molecular_weight": round(mol_wt, 4)
        }
    except Exception as e:
        print(f"Error logic: {e}")
        return None

if __name__ == '__main__':
    smiles = 'CCO'
    result = level_function(smiles)
    if result:
        print(f'Output: {result}')
