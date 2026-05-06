from rdkit import Chem

def level_function(mol):
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        atoms = mol.GetAtoms()
        oxygen_count = sum(1 for atom in atoms if atom.GetSymbol() == "O")
        return oxygen_count
    except Exception as e:
        print(e)
        return None
