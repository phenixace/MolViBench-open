from rdkit import Chem

def level_function(mol):
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return set([atom.GetSymbol() for atom in mol.GetAtoms()])
    except Exception as e:
        print(e)
        return None
