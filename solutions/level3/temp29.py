from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol_smiles):



    try:
        molecule = Chem.MolFromSmiles(mol_smiles)
        if molecule is None:
            return None





        pattern = '[C;H1,H2,H3:1]-[C:2]-[F,Cl,Br,I:3]>>[C:1]=[C:2]'
        rxn = AllChem.ReactionFromSmarts(pattern)


        products = rxn.RunReactants((molecule,))

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
    smiles = 'CC(Br)CC'
    print(f'Output: {smiles}')
    print(f'Output: {level_function(smiles)}')
