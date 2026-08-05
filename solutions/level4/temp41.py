from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol_smi):



    try:
        mol_obj = Chem.MolFromSmiles(mol_smi)
        if mol_obj is None:
            return None



        pattern = Chem.MolFromSmarts('[CX3H1](=O)')
        has_aldehyde = mol_obj.HasSubstructMatch(pattern)

        if not has_aldehyde:
            return None





        rxn_smarts = '[*:2][CX3H1:1](=O)>>[*:2][CX4H2:1][OH]'
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)

        products = rxn.RunReactants((mol_obj,))

        if not products:

            rxn_fallback = AllChem.ReactionFromSmarts('[CH1:1](=O)>>[CH2:1]O')
            products = rxn_fallback.RunReactants((mol_obj,))

        if not products:
            return {"has_aldehyde": True, "product": "Reduction failed", "logp": None}


        product = products[0][0]


        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)



        logp = rdMolDescriptors.CalcCrippenDescriptors(product)[0]

        return {
            "has_aldehyde": True,
            "product": product_smiles,
            "logp": round(logp, 4)
        }
    except Exception as e:
        print(f"Error logic: {e}")
        return None

if __name__ == '__main__':
    smiles = 'CC=O'
    result = level_function(smiles)
    if result:
        print(f'Output: {result}')
