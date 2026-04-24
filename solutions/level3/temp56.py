from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold


def level_function(mol, new_scaffold_smiles):
    """进行 scaffold hopping：保留分子侧链，替换核心骨架为给定新骨架。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        new_scaffold = Chem.MolFromSmiles(new_scaffold_smiles)
        if mol_obj is None or new_scaffold is None:
            return None

        # Extract original scaffold
        original_scaffold = MurckoScaffold.GetScaffoldForMol(mol_obj)
        original_scaffold_smi = Chem.MolToSmiles(original_scaffold)

        # Get side chains by removing scaffold
        # Use ReplaceSubstructs to replace the scaffold
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


if __name__ == "__main__":
    mol_smi = "c1ccc(NC(=O)C)cc1"
    new_scaf = "c1ccncc1"
    result = level_function(mol_smi, new_scaf)
    if result:
        print(f"原始骨架: {result['original_scaffold']}")
        print(f"新骨架: {result['new_scaffold']}")
        print(f"替换结果: {result['hopped_molecules']}")
