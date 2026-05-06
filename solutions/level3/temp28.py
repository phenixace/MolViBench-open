from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol1, mol2):
    try:
        carbonyl = Chem.MolFromSmiles(mol1)
        nucleophile = Chem.MolFromSmiles(mol2)
        if carbonyl is None or nucleophile is None:
            return None
        reaction_smarts = '[C:1](=[O:2]).[Nu:3]>>[C:1]([O:2])[Nu:3]'
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        products = rxn.RunReactants((carbonyl, nucleophile))
        if not products:
            reaction_smarts_alt = '[C:1](=[O:2])([C:3])[C:4].[N:5]>>[C:1]([OH:2])([C:3])([C:4])[N:5]'
            rxn = AllChem.ReactionFromSmarts(reaction_smarts_alt)
            products = rxn.RunReactants((carbonyl, nucleophile))
        if not products:
            reaction_smarts_simple = '[C:1]=[O:2].[N:3]>>[C:1]([OH:2])[N:3]'
            rxn = AllChem.ReactionFromSmarts(reaction_smarts_simple)
            products = rxn.RunReactants((carbonyl, nucleophile))
        if not products:
            return None
        result_smiles = []
        for product_set in products:
            for product in product_set:
                try:
                    Chem.SanitizeMol(product)
                    smi = Chem.MolToSmiles(product)
                    if smi not in result_smiles:
                        result_smiles.append(smi)
                except Exception:
                    continue
        return result_smiles if result_smiles else None
    except Exception as e:
        print(e)
        return None
