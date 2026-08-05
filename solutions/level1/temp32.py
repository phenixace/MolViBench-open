from rdkit import Chem

def level_function(mol):



    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        atom_features = []
        for atom in mol.GetAtoms():
            features = [
                atom.GetAtomicNum(),
                atom.GetDegree(),
                atom.GetFormalCharge(),
                atom.GetHybridization().real,
                int(atom.GetIsAromatic())
            ]
            atom_features.append(features)
        return atom_features
    except Exception as e:
        print(e)
        return None

if __name__ == '__main__':
    smiles = 'CC[C@H](F)C(=O)O'
    print(f'Output: {level_function(smiles)}')
