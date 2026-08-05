from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(azide_smiles, alkyne_smiles):
    try:
        azide = Chem.MolFromSmiles(azide_smiles)
        alkyne = Chem.MolFromSmiles(alkyne_smiles)
        if azide is None or alkyne is None:
            return None










        rxn_smarts = '[C:1]#[C:2].[N:3]=[N+:4]=[N-:5]>>[N:3]1-[N:4]=[N:5]-[C:1]=[C:2]1'

        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        products = rxn.RunReactants((alkyne, azide))

        if not products:
            return None

        results = set()
        for prod_set in products:
            for prod in prod_set:
                try:

                    for atom in prod.GetAtoms():
                        atom.SetFormalCharge(0)


                    Chem.SanitizeMol(prod)
                    results.add(Chem.MolToSmiles(prod))
                except Exception:
                    continue

        return sorted(list(results)) if results else None
    except Exception as e:
        print(f"Error: {e}")
        return None

if __name__ == '__main__':
    azide = 'c1ccc(N=[N+]=[N-])cc1'
    alkyne = 'C#Cc1ccccc1'
    result = level_function(azide, alkyne)
    print(f'Output: {result}')
