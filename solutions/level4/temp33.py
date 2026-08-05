from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol_smi):



    try:
        mol_obj = Chem.MolFromSmiles(mol_smi)
        if mol_obj is None:
            return None



        pattern = Chem.MolFromSmarts('[CX4][OX2H]')
        has_alcohol = mol_obj.HasSubstructMatch(pattern)

        if not has_alcohol:
            return None





        rxn_smarts = '[*:1][CX4:2][CX4:3]([OX2H])[*:4]>>[*:1][CX3:2]=[CX3:3][*:4]'
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)

        products = rxn.RunReactants((mol_obj,))
        if not products:

            rxn_simple = AllChem.ReactionFromSmarts('[CX4:1][CX4:2][OH]>>[CX3:1]=[CX3:2]')
            products = rxn_simple.RunReactants((mol_obj,))

        if not products:
            return {"has_alcohol": True, "product": "Dehydration failed", "molecular_weight": None}


        product = products[0][0]


        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)


        mol_wt = rdMolDescriptors.CalcExactMolWt(product)

        return {
            "has_alcohol": True,
            "product": product_smiles,
            "molecular_weight": round(mol_wt, 4)
        }
    except Exception as e:
        print(f"Error logic: {e}")
        return None

if __name__ == '__main__':
    smiles = 'CCCO'
    result = level_function(smiles)
    if result:
        print(f'Output: {result}')
