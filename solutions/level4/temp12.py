from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, RWMol

def level_function(mol):
    """给定分子 → 判断是否含杂原子 → 若有 → 替换为氧 → 计算 TPSA。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含杂原子（非 C、H 的原子）
        heteroatom_idx = None
        for atom in mol_obj.GetAtoms():
            if atom.GetAtomicNum() not in (6, 1):
                heteroatom_idx = atom.GetIdx()
                break

        has_heteroatom = heteroatom_idx is not None

        if not has_heteroatom:
            return None

        # Step 2: 将第一个杂原子替换为氧
        rw_mol = RWMol(mol_obj)
        rw_mol.GetAtomWithIdx(heteroatom_idx).SetAtomicNum(8)

        product = rw_mol.GetMol()
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        # Step 3: 计算 TPSA
        tpsa = rdMolDescriptors.CalcTPSA(product)

        return {
            "has_heteroatom": has_heteroatom,
            "product": product_smiles,
            "tpsa": round(tpsa, 4)
        }
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CCN"
    print(f"result: {level_function(smiles)}")
