from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

def level_function(mol_smi):



    try:
        mol_obj = Chem.MolFromSmiles(mol_smi)
        if mol_obj is None:
            return None


        pattern = Chem.MolFromSmarts('[OX2H]')
        has_hydroxyl = mol_obj.HasSubstructMatch(pattern)

        if not has_hydroxyl:
            return None




        rxn_smarts = '[*:1][OH:2]>>[*:1]Cl'
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)

        products = rxn.RunReactants((mol_obj,))

        if not products:

            rxn_fallback = AllChem.ReactionFromSmarts('[C:1][OH]>>[C:1]Cl')
            products = rxn_fallback.RunReactants((mol_obj,))

        if not products:
            return {"has_hydroxyl": True, "product": "Halogenation failed", "qed": None}


        product = products[0][0]


        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)


        qed = Descriptors.qed(product)

        return {
            "has_hydroxyl": True,
            "product": product_smiles,
            "qed": round(qed, 4)
        }
    except Exception as e:
        print(f"Error logic: {e}")
        return None

if __name__ == '__main__':
    smiles = 'c1ccc(O)cc1'
    result = level_function(smiles)
    if result:
        print(f'Output: {result}')
