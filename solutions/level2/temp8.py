from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        rxn = AllChem.ReactionFromSmarts('[cH:1]>>[c:1]C')
        first_products = rxn.RunReactants((mol_obj,))
        first_unique = set()
        first_mols = []
        for product_tuple in first_products:
            for product in product_tuple:
                Chem.SanitizeMol(product)
                smi = Chem.MolToSmiles(product)
                if smi not in first_unique:
                    first_unique.add(smi)
                    first_mols.append(product)
        unique_smiles = set()
        for m in first_mols:
            second_products = rxn.RunReactants((m,))
            for product_tuple in second_products:
                for product in product_tuple:
                    Chem.SanitizeMol(product)
                    smi = Chem.MolToSmiles(product)
                    unique_smiles.add(smi)
        return sorted(list(unique_smiles))
    except Exception as e:
        print(e)
        return None
