from rdkit import Chem

def level_function(mol):
    """
    判断分子是否含有五元环。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if len(ring) == 5:
                return True
        return False
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CCOc1ccccc1"
    print(f"是否含有五元环: {level_function(smiles)}")