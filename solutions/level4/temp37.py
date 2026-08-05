from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol_smi):
    try:
        mol_obj = Chem.MolFromSmiles(mol_smi)
        if mol_obj is None:
            return None


        pattern = Chem.MolFromSmarts('[#6][CX3](=O)[#6]')
        if not mol_obj.HasSubstructMatch(pattern):
            return None



        mol_with_hs = Chem.AddHs(mol_obj)



        rxn_smarts = '[C:1](-[H:4])-[C:2]=[O:3]>>[C:1]=[C:2]-[O:3]-[H:4]'
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)

        products = rxn.RunReactants((mol_with_hs,))

        if not products:
            return {"has_ketone": True, "product": "Enolization failed", "formula": None}


        product_with_hs = products[0][0]
        product = Chem.RemoveHs(product_with_hs)


        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)
        formula = rdMolDescriptors.CalcMolFormula(product)

        return {
            "has_ketone": True,
            "product": product_smiles,
            "formula": formula
        }

    except Exception as e:
        print(f"Error logic: {e}")
        return None

if __name__ == '__main__':
    smiles = 'CC(=O)CC'
    result = level_function(smiles)
    if result:
        print(f'Output: {result}')
