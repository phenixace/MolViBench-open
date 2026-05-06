from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        remover = rdMolStandardize.LargestFragmentChooser()
        mol_obj = remover.choose(mol_obj)

        te = rdMolStandardize.TautomerEnumerator()
        mol_obj = te.Canonicalize(mol_obj)

        uncharger = rdMolStandardize.Uncharger()
        mol_obj = uncharger.uncharge(mol_obj)

        Chem.SanitizeMol(mol_obj)
        return Chem.MolToSmiles(mol_obj)
    except Exception as e:
        print(e)
        return None
