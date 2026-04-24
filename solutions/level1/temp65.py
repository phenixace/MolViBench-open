from rdkit import Chem
from rdkit.Chem import Crippen


def level_function(mol):
    """计算分子的 Crippen 分子折射率 (MR) 的原子贡献分解。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        contribs = Crippen._GetAtomContribs(mol_obj)
        # contribs is list of (logP_contrib, MR_contrib)
        atom_mr = {i: round(c[1], 4) for i, c in enumerate(contribs)}
        total_mr = round(sum(c[1] for c in contribs), 4)
        return {"atom_contributions": atom_mr, "total_MR": total_mr}
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "CCO"
    print(f"Crippen MR 分解: {level_function(smiles)}")
