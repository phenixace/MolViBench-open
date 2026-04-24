from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):
    """给定分子 → 判断是否含芳香环 → 若有 → 生成构象 → 优化能量。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含芳香环
        num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol_obj)
        has_aromatic_ring = num_aromatic_rings > 0

        if not has_aromatic_ring:
            return None

        # Step 2: 生成3D构象并优化能量
        mol_h = Chem.AddHs(mol_obj)
        embed_result = AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())
        if embed_result == -1:
            return None

        optimize_result = AllChem.MMFFOptimizeMolecule(mol_h)

        # Step 3: 计算优化后能量
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

if __name__ == "__main__":
    smiles = "c1ccccc1"
    print(f"result: {level_function(smiles)}")
