from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def level_function(mol):
    """给定分子 → 判断是否含杂环 → 若有 → 打开环 → 再计算分子量。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含杂环
        ring_info = mol_obj.GetRingInfo()
        has_heterocycle = False
        target_bond = None

        for ring in ring_info.BondRings():
            atom_rings = ring_info.AtomRings()
            for atom_ring in atom_rings:
                if any(mol_obj.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in atom_ring):
                    has_heterocycle = True
                    # 找到杂原子参与的键用于断裂
                    for bond_idx in ring:
                        bond = mol_obj.GetBondWithIdx(bond_idx)
                        a1 = bond.GetBeginAtom()
                        a2 = bond.GetEndAtom()
                        if a1.GetAtomicNum() != 6 or a2.GetAtomicNum() != 6:
                            target_bond = bond_idx
                            break
                    break
            if has_heterocycle:
                break

        if not has_heterocycle:
            return None

        # Step 2: 打开环 (断裂含杂原子的键)
        rw = Chem.RWMol(mol_obj)
        if target_bond is not None:
            bond = rw.GetBondWithIdx(target_bond)
            rw.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        else:
            # 如果没找到杂原子键, 断裂环中第一个键
            bonds = list(ring_info.BondRings())[0]
            bond = rw.GetBondWithIdx(bonds[0])
            rw.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

        try:
            Chem.SanitizeMol(rw)
        except Exception:
            pass
        product_smiles = Chem.MolToSmiles(rw)

        # Step 3: 计算分子量
        product_mol = Chem.MolFromSmiles(product_smiles)
        if product_mol is None:
            return None
        mol_wt = rdMolDescriptors.CalcExactMolWt(product_mol)

        return {
            "has_heterocycle": has_heterocycle,
            "product": product_smiles,
            "molecular_weight": round(mol_wt, 4)
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccncc1"  # 吡啶
    print(f"result: {level_function(smiles)}")
