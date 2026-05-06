from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol_obj)
        has_aromatic_ring = num_aromatic_rings > 0

        if not has_aromatic_ring:
            return None

        mol_h = Chem.AddHs(mol_obj)
        embed_result = AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())
        if embed_result == -1:
            return None

        optimize_result = AllChem.MMFFOptimizeMolecule(mol_h)

        ff = AllChem.MMFFGetMoleculeForceField(mol_h, AllChem.MMFFGetMoleculeProperties(mol_h))
        if ff is None:
            return None
        energy = ff.CalcEnergy()

        product_smiles = Chem.MolToSmiles(mol_obj)

        return {
            "has_aromatic_ring": has_aromatic_ring,
            "product": product_smiles,
            "energy": round(energy, 4)
        }
    except Exception as e:
        print(e)
        return None
