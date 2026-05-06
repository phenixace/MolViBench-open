from rdkit import Chem

def level_function(mol):
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        atom_features = []
        for atom in mol.GetAtoms():
            features = [
            ]
            atom_features.append(features)
        return atom_features
    except Exception as e:
        print(e)
        return None
