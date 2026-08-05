from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def level_function(mol_smi):
    try:
        mol = Chem.MolFromSmiles(mol_smi)
        if mol is None: return None



        hetero_in_ring = mol.HasSubstructMatch(Chem.MolFromSmarts("[!#6;r]"))

        if not hetero_in_ring:
            return None


        rw = Chem.RWMol(mol)
        ring_info = rw.GetRingInfo()


        for atom_ring in ring_info.AtomRings():
            has_hetero = any(rw.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in atom_ring)
            if has_hetero:

                found_bond = False
                for i in range(len(atom_ring)):
                    a1, a2 = atom_ring[i], atom_ring[(i+1)%len(atom_ring)]
                    bond = rw.GetBondBetweenAtoms(a1, a2)
                    if rw.GetAtomWithIdx(a1).GetAtomicNum() != 6 or rw.GetAtomWithIdx(a2).GetAtomicNum() != 6:
                        rw.RemoveBond(a1, a2)
                        found_bond = True
                        break
                if found_bond: break


        for atom in rw.GetAtoms():
            atom.SetIsAromatic(False)
            atom.SetNumExplicitHs(0)

        for bond in rw.GetBonds():
            bond.SetIsAromatic(False)
            if bond.GetBondType() == Chem.rdchem.BondType.AROMATIC:
                bond.SetBondType(Chem.rdchem.BondType.SINGLE)


        Chem.SanitizeMol(rw)

        final_mol = Chem.AddHs(rw)
        final_mol = Chem.RemoveHs(final_mol)

        product_smiles = Chem.MolToSmiles(final_mol)
        mol_wt = rdMolDescriptors.CalcExactMolWt(final_mol)

        return {
            "has_heterocycle": True,
            "product": product_smiles,
            "molecular_weight": round(mol_wt, 4)
        }
    except Exception as e:
        print(f"Error: {e}")
        return None

if __name__ == '__main__':
    print(f"Output: {level_function('c1ccncc1')}")
