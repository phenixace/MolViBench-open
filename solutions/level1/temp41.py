from rdkit import Chem

def level_function(mol):
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        atoms = mol.GetAtoms()
        for atom in atoms:
            if atom.GetSymbol() == "P":
                return True
        return False
    except Exception as e:
        print(e)
        return None
