from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


def level_function(mol, num_confs=50):
    """生成大环分子的多个 3D 构象并评估环张力。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Check if it has a macrocycle (ring >= 12)
        ring_info = mol_obj.GetRingInfo()
        max_ring_size = max((len(r) for r in ring_info.AtomRings()), default=0)

        mol_h = Chem.AddHs(mol_obj)

        # Generate conformers
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.numThreads = 1
        params.useRandomCoords = True if max_ring_size >= 12 else False

        cids = AllChem.EmbedMultipleConfs(mol_h, numConfs=num_confs, params=params)
        if not cids:
            return None

        # Optimize with MMFF and collect energies
        energies = []
        for cid in cids:
            try:
                result = AllChem.MMFFOptimizeMolecule(mol_h, confId=cid, maxIters=500)
                ff = AllChem.MMFFGetMoleculeForceField(mol_h, AllChem.MMFFGetMoleculeProperties(mol_h), confId=cid)
                if ff:
                    energy = ff.CalcEnergy()
                    energies.append({"conf_id": cid, "energy_kcal": round(energy, 2)})
            except Exception:
                continue

        if not energies:
            return None

        energies.sort(key=lambda x: x["energy_kcal"])
        min_e = energies[0]["energy_kcal"]
        max_e = energies[-1]["energy_kcal"]

        return {
            "max_ring_size": max_ring_size,
            "is_macrocycle": max_ring_size >= 12,
            "num_conformers": len(energies),
            "min_energy": min_e,
            "max_energy": max_e,
            "energy_range": round(max_e - min_e, 2),
            "top3_conformers": energies[:3]
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    # Cyclododecane - a 12-membered macrocycle
    smiles = "C1CCCCCCCCCCC1"
    result = level_function(smiles, num_confs=20)
    if result:
        print(f"Macrocycle: {result['is_macrocycle']}, Ring size: {result['max_ring_size']}")
        print(f"Energy range: {result['energy_range']} kcal/mol")
