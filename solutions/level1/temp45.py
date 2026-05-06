from rdkit import Chem

def level_function(mol):
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        atoms = mol.GetAtoms()
        nitrogen_count = sum(1 for atom in atoms if atom.GetSymbol() == "N")
        return nitrogen_count
    except Exception as e:
        print(e)
        return None
