from rdkit import Chem


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        ring_info = mol_obj.GetRingInfo()
        for ring in ring_info.AtomRings():
            if len(ring) >= 12:
                return True
        return False
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'C1CCCCCCCCCCCCC1'
    print(f'Output: {level_function(smiles)}')
