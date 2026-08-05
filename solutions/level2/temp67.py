from rdkit import Chem


def level_function(smarts_pattern, library_smiles):

    try:
        pattern = Chem.MolFromSmarts(smarts_pattern)
        if pattern is None:
            return None

        results = []
        for smi in library_smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            matches = mol.GetSubstructMatches(pattern)
            if matches:
                results.append({
                    "smiles": Chem.MolToSmiles(mol),
                    "matches": [list(m) for m in matches]
                })

        return results
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    pattern = '[OH]'
    library = ['c1ccccc1', 'c1ccc(O)cc1', 'CCO', 'CC(=O)O', 'CCCC']
    result = level_function(pattern, library)
    for r in result:
        print(f"Output: {r['smiles']}{r['matches']}")
