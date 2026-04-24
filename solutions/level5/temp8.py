from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions


def level_function(mol):
    """给定一个小分子，生成其立体异构体并保留所有可行构象。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Enumerate stereoisomers
        opts = StereoEnumerationOptions(tryEmbedding=True, unique=True)
        stereoisomers = list(EnumerateStereoisomers(mol_obj, options=opts))

        results = {}
        for iso in stereoisomers:
            iso_smi = Chem.MolToSmiles(iso)

            # Generate 3D conformers for each stereoisomer
            iso_3d = Chem.AddHs(iso)
            try:
                params = AllChem.EmbedMultipleConfs(
                    iso_3d,
                    numConfs=10,
                    randomSeed=42,
                    pruneRmsThresh=0.5,
                    maxAttempts=100
                )
                num_confs = iso_3d.GetNumConformers()

                # Optimize conformers
                if num_confs > 0:
                    try:
                        AllChem.MMFFOptimizeMoleculeConfs(iso_3d)
                    except Exception:
                        pass
                    num_confs = iso_3d.GetNumConformers()

                results[iso_smi] = num_confs
            except Exception:
                results[iso_smi] = 0

        return results
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    result = level_function("CC(O)C(F)Cl")
    print(f"result: {result}")
