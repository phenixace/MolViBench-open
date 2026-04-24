from rdkit import Chem


def level_function(mol):
    """输出分子中所有双键的 E/Z 构型。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        Chem.AssignStereochemistry(mol_obj, cleanIt=True, force=True)
        ez_bonds = []
        for bond in mol_obj.GetBonds():
            stereo = bond.GetStereo()
            if stereo != Chem.BondStereo.STEREONONE:
                begin_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()
                stereo_str = str(stereo).split('.')[-1]
                ez_bonds.append((begin_idx, end_idx, stereo_str))
        return ez_bonds
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = r"C/C=C\C"
    print(f"E/Z 构型: {level_function(smiles)}")
