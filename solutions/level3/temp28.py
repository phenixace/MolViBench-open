from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(carbonyl_smi, nucleophile_smi):



    try:
        carbonyl = Chem.MolFromSmiles(carbonyl_smi)
        nucleophile = Chem.MolFromSmiles(nucleophile_smi)
        if carbonyl is None or nucleophile is None:
            return None





        rxn_smarts = '[C:1]=[O:2].[N,O,S;H1,H2,H3:3]>>[C:1](-[O;H1:2])-[*:3]'

        rxn = AllChem.ReactionFromSmarts(rxn_smarts)


        products = rxn.RunReactants((carbonyl, nucleophile))

        result_smiles = set()
        if products:
            for product_set in products:
                for mol in product_set:
                    try:

                        Chem.SanitizeMol(mol)
                        smi = Chem.MolToSmiles(mol)
                        result_smiles.add(smi)
                    except:
                        continue

        return list(result_smiles) if result_smiles else None
    except Exception as e:
        print(f"Error: {e}")
        return None

if __name__ == '__main__':
    print(f"Output: {level_function('CC=O', 'CN')}")
