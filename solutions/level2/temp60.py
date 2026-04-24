from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


def level_function(mol):
    """对分子进行标准化处理（去盐、标准化互变异构体、标准化电荷）。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: Remove salts - keep largest fragment
        remover = rdMolStandardize.LargestFragmentChooser()
        mol_obj = remover.choose(mol_obj)

        # Step 2: Standardize tautomer
        te = rdMolStandardize.TautomerEnumerator()
        mol_obj = te.Canonicalize(mol_obj)

        # Step 3: Uncharge (neutralize)
        uncharger = rdMolStandardize.Uncharger()
        mol_obj = uncharger.uncharge(mol_obj)

        # Step 4: Cleanup
        Chem.SanitizeMol(mol_obj)
        return Chem.MolToSmiles(mol_obj)
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "[Na+].[O-]c1ccccc1"
    print(f"标准化结果: {level_function(smiles)}")
