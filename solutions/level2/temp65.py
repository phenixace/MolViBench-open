from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(core_smiles, rgroup_lists):
    """给定核心骨架和多组 R-group 的 SMILES 列表，枚举所有组合产物并返回 canonical SMILES 列表。

    参数:
        core_smiles: 骨架 SMILES，R-group 位置用 [*:1], [*:2] 等标记
        rgroup_lists: 字典，key 为 R-group 编号 (1,2,...), value 为 SMILES 列表
    """
    try:
        from itertools import product

        core = Chem.MolFromSmiles(core_smiles)
        if core is None:
            return None

        # Get dummy atom indices in core
        dummy_info = {}
        for atom in core.GetAtoms():
            if atom.GetAtomicNum() == 0:  # dummy atom
                map_num = atom.GetAtomMapNum()
                if map_num > 0:
                    dummy_info[map_num] = atom.GetIdx()

        if not dummy_info:
            return None

        # Prepare R-group combinations
        rgroup_keys = sorted(rgroup_lists.keys())
        rgroup_values = [rgroup_lists[k] for k in rgroup_keys]

        results = set()
        for combo in product(*rgroup_values):
            try:
                combined = Chem.RWMol(core)
                # Process each R-group attachment
                for rg_idx, rg_smi in zip(rgroup_keys, combo):
                    rg_mol = Chem.MolFromSmiles(rg_smi)
                    if rg_mol is None:
                        continue

                    # Use reaction-based attachment
                    rxn_smarts = f"[*:{rg_idx}][*:99].[*:99]{rg_smi}>>[*:{rg_idx}]{rg_smi}"

                # Alternative: use AllChem.ReplaceSubstructs
                mol = Chem.RWMol(core)
                for rg_idx, rg_smi in zip(rgroup_keys, combo):
                    rg_mol = Chem.MolFromSmiles(rg_smi)
                    if rg_mol is None:
                        continue
                    # Find dummy with this map number
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


if __name__ == "__main__":
    core = "[*:1]c1ccc([*:2])cc1"
    rgroups = {1: ["C", "CC", "O"], 2: ["F", "Cl", "N"]}
    result = level_function(core, rgroups)
    print(f"组合产物数: {len(result) if result else 0}")
    if result:
        for smi in result[:5]:
            print(f"  {smi}")
