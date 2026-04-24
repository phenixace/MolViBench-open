from rdkit import Chem


def level_function(mol):
    """判断分子是否为大环分子（含有 ≥12 元环）。"""
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


if __name__ == "__main__":
    smiles = "C1CCCCCCCCCCCCC1"  # 14-membered ring
    print(f"大环分子: {level_function(smiles)}")
