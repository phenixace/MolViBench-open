from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None


        opts = StereoEnumerationOptions(tryEmbedding=True, unique=True)
        stereoisomers = list(EnumerateStereoisomers(mol_obj, options=opts))

        results = {}
        for iso in stereoisomers:
            iso_smi = Chem.MolToSmiles(iso)


            iso_3d = Chem.AddHs(iso)
            try:
                if iso_3d.GetNumHeavyAtoms() > 50:
                    results[iso_smi] = 0
                    continue
                params = AllChem.EmbedMultipleConfs(
                    iso_3d,
                    numConfs=10,
                    randomSeed=42,
                    pruneRmsThresh=0.5,
                    maxAttempts=5
                )
                num_confs = iso_3d.GetNumConformers()


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


if __name__ == '__main__':
    result = level_function('CC(O)C(F)Cl')
    print(f'Output: {result}')
