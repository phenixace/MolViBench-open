from rdkit import Chem


def level_function(mol, filename="output.sdf"):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        writer = Chem.SDWriter(filename)
        writer.write(mol_obj)
        writer.close()
        return filename
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'CCO'
    result = level_function(smiles, 'output.sdf')
    print(f'Output: {result}')
