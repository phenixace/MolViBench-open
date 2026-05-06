from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(core_smiles, rgroup_lists):
    try:
        from itertools import product

        core = Chem.MolFromSmiles(core_smiles)
        if core is None:
            return None

        dummy_info = {}
        for atom in core.GetAtoms():
                map_num = atom.GetAtomMapNum()
                if map_num > 0:
                    dummy_info[map_num] = atom.GetIdx()

        if not dummy_info:
            return None

        rgroup_keys = sorted(rgroup_lists.keys())
        rgroup_values = [rgroup_lists[k] for k in rgroup_keys]

        results = set()
        for combo in product(*rgroup_values):
            try:
                combined = Chem.RWMol(core)
                for rg_idx, rg_smi in zip(rgroup_keys, combo):
                    rg_mol = Chem.MolFromSmiles(rg_smi)
                    if rg_mol is None:
                        continue

                    rxn_smarts = f"[*:{rg_idx}][*:99].[*:99]{rg_smi}>>[*:{rg_idx}]{rg_smi}"

                mol = Chem.RWMol(core)
                for rg_idx, rg_smi in zip(rgroup_keys, combo):
                    rg_mol = Chem.MolFromSmiles(rg_smi)
                    if rg_mol is None:
                        continue
                    dummy_pattern = Chem.MolFromSmarts(f"[#0:{rg_idx}]")
                    replaced = AllChem.ReplaceSubstructs(mol, dummy_pattern, rg_mol,
                                                         replaceAll=True)
                    if replaced:
                        mol = Chem.RWMol(replaced[0])

                try:
                    Chem.SanitizeMol(mol)
                    smi = Chem.MolToSmiles(mol)
                    if smi:
                        results.add(smi)
                except Exception:
                    pass
            except Exception:
                continue

        return sorted(list(results))
    except Exception as e:
        print(e)
        return None
