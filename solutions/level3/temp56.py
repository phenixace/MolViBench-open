from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def level_function(mol, new_scaffold_smiles):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        new_scaffold = Chem.MolFromSmiles(new_scaffold_smiles)
        if mol_obj is None or new_scaffold is None:
            return None

        original_scaffold = MurckoScaffold.GetScaffoldForMol(mol_obj)
        original_scaffold_smi = Chem.MolToSmiles(original_scaffold)

        replaced = AllChem.ReplaceSubstructs(mol_obj, original_scaffold, new_scaffold)
        if not replaced:
            return None

        results = []
        for prod in replaced:
            try:
                Chem.SanitizeMol(prod)
                smi = Chem.MolToSmiles(prod)
                if smi:
                    results.append(smi)
            except Exception:
                continue

        if not results:
            return None

        return {
            "original_scaffold": original_scaffold_smi,
            "new_scaffold": Chem.MolToSmiles(new_scaffold),
            "hopped_molecules": list(set(results))
        }
    except Exception as e:
        print(e)
        return None
