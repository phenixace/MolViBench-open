from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol):
    try:
        molecule = Chem.MolFromSmiles(mol)
        if molecule is None:
            return None

        molecule = Chem.AddHs(molecule)
        AllChem.EmbedMolecule(molecule, randomSeed=42)

        num_confs = 20
        cids = AllChem.EmbedMultipleConfs(
            molecule,
            numConfs=num_confs,
            randomSeed=42,
            pruneRmsThresh=0.5,
            useExpTorsionAnglePrefs=True,
            useBasicKnowledge=True,
        )

        if not cids:
            return None

        results = AllChem.MMFFOptimizeMoleculeConfs(molecule)

        conformers = []
        for i, cid in enumerate(cids):
            energy = results[i][1] if results[i][0] == 0 else None
            molblock = Chem.MolToMolBlock(molecule, confId=cid)
            conformers.append({
                "conf_id": cid,
                "energy": energy,
                "molblock": molblock,
            })

        conformers.sort(key=lambda x: x["energy"] if x["energy"] is not None else float("inf"))
        return conformers
    except Exception as e:
        print(e)
        return None
