from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        ring_info = mol_obj.GetRingInfo()
        has_heterocycle = False
        target_bond = None

        for ring in ring_info.BondRings():
            atom_rings = ring_info.AtomRings()
            for atom_ring in atom_rings:
                if any(mol_obj.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in atom_ring):
                    has_heterocycle = True
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

        rw = Chem.RWMol(mol_obj)
        if target_bond is not None:
            bond = rw.GetBondWithIdx(target_bond)
            rw.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        else:
            bonds = list(ring_info.BondRings())[0]
            bond = rw.GetBondWithIdx(bonds[0])
            rw.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

        try:
            Chem.SanitizeMol(rw)
        except Exception:
            pass
        product_smiles = Chem.MolToSmiles(rw)

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
