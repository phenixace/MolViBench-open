from rdkit import Chem


def level_function(mol):
    """获取分子中每个原子的杂化类型（sp, sp2, sp3）。"""
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


if __name__ == "__main__":
    smiles = "C=CC#N"
    print(f"杂化类型: {level_function(smiles)}")
