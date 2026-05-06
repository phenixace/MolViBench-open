from rdkit import Chem

def level_function(mol):
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in ring):
                return True
        return False
    except Exception as e:
        print(e)
        return None
