from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, RWMol

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        ring_info = mol_obj.GetRingInfo()
        has_ring = ring_info.NumRings() > 0

        if not has_ring:
            return None

        rw_mol = Chem.RWMol(mol_obj)
        bond_rings = ring_info.BondRings()
        if not bond_rings:
            return None

        first_ring_bonds = bond_rings[0]
        bond_idx = first_ring_bonds[0]
        bond = rw_mol.GetBondWithIdx(bond_idx)
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        rw_mol.RemoveBond(begin_idx, end_idx)

        try:
            Chem.SanitizeMol(rw_mol)
        except Exception:
            pass

        product_smiles = Chem.MolToSmiles(rw_mol)

        product = Chem.MolFromSmiles(product_smiles)
        if product is None:
            mol_wt = rdMolDescriptors.CalcExactMolWt(rw_mol)
        else:
            mol_wt = rdMolDescriptors.CalcExactMolWt(product)

        return {
            "has_ring": has_ring,
            "product": product_smiles,
            "molecular_weight": round(mol_wt, 4)
        }
    except Exception as e:
        print(e)
        return None
