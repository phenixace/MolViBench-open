from rdkit import Chem


def level_function(mol):
    """给定分子 → 判断是否含芳香氮 → 若有 → 质子化 → 计算电荷数。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含芳香氮
        pattern = Chem.MolFromSmarts("n")
        has_aromatic_n = mol_obj.HasSubstructMatch(pattern)

        if not has_aromatic_n:
            return None

        # Step 2: 质子化芳香氮
        rw = Chem.RWMol(mol_obj)
        for atom in rw.GetAtoms():
            if atom.GetIsAromatic() and atom.GetAtomicNum() == 7:
                atom.SetFormalCharge(1)
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
                break  # 只质子化第一个芳香氮

        try:
            Chem.SanitizeMol(rw)
        except Exception:
            pass
        product_smiles = Chem.MolToSmiles(rw)

        # Step 3: 计算电荷数
        total_charge = sum(atom.GetFormalCharge() for atom in rw.GetAtoms())

        return {
            "has_aromatic_n": has_aromatic_n,
            "product": product_smiles,
            "total_charge": total_charge
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccncc1"  # 吡啶
    print(f"result: {level_function(smiles)}")
