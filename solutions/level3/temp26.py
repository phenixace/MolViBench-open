from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol):

    try:
        molecule = Chem.MolFromSmiles(mol)
        if molecule is None:
            return None
        reaction_patterns = [
            '[CH3:1]>>[CH2:1]Cl',
            '[CH2:1]>>[CH:1]Cl',
            '[CH:1]>>[C:1]Cl',
        ]
        result_smiles = []
        for pattern in reaction_patterns:
            rxn = AllChem.ReactionFromSmarts(pattern)
            products = rxn.RunReactants((molecule,))
            if products:
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


if __name__ == '__main__':
    smiles = 'CCC'
    print(f'Output: {level_function(smiles)}')
