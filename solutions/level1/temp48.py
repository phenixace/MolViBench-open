from rdkit import Chem

def level_function(mol):
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        atoms = mol.GetAtoms()
        halogens = {"F", "Cl", "Br", "I", "At"}
        return sum(1 for atom in atoms if atom.GetSymbol() in halogens)
    except Exception as e:
        print(e)
        return None
