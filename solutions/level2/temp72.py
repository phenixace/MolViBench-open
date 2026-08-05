from rdkit import Chem


def level_function(mol, smarts_pattern):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        pattern = Chem.MolFromSmarts(smarts_pattern)
        if pattern is None:
            return None

        matches = mol_obj.GetSubstructMatches(pattern)
        return [list(m) for m in matches]
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccc(O)c(O)c1'
    smarts = '[OH]'
    result = level_function(smiles, smarts)
    print(f'Output: {result}')
