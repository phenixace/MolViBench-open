from rdkit import Chem

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        hybridizations = {}
        for atom in mol_obj.GetAtoms():
            hyb = atom.GetHybridization()
            hybridizations[atom.GetIdx()] = str(hyb).split('.')[-1]
        return hybridizations
    except Exception as e:
        print(e)
        return None
