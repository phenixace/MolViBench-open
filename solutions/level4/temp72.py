from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, rdMolDescriptors
from itertools import product


def level_function(core_smiles, rgroup_dict, top_k=5):
    """给定核心骨架和 3 组 R-group → 枚举所有组合产物 → Lipinski 过滤 → 按 QED 降序排列 → 输出 Top-5 分子及其 QED 值。"""
    try:
        core = Chem.MolFromSmiles(core_smiles)
        if core is None:
            return None

        rg_keys = sorted(rgroup_dict.keys())
        rg_values = [rgroup_dict[k] for k in rg_keys]

        all_products = set()
        for combo in product(*rg_values):
            mol = Chem.RWMol(core)
            for rg_idx, rg_smi in zip(rg_keys, combo):
                rg_mol = Chem.MolFromSmiles(rg_smi)
                if rg_mol is None:
                    continue
                dummy_pattern = Chem.MolFromSmarts(f"[#0:{rg_idx}]")
                replaced = AllChem.ReplaceSubstructs(mol, dummy_pattern, rg_mol, replaceAll=True)
                if replaced:
                    mol = Chem.RWMol(replaced[0])
            try:
                Chem.SanitizeMol(mol)
                all_products.add(Chem.MolToSmiles(mol))
            except Exception:
                pass

        # Lipinski filter + QED
        scored = []
        for smi in all_products:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            hbd = rdMolDescriptors.CalcNumHBD(mol)
            hba = rdMolDescriptors.CalcNumHBA(mol)
            if mw < 500 and logp < 5 and hbd <= 5 and hba <= 10:
                qed = Descriptors.qed(mol)
                scored.append({"smiles": smi, "QED": round(qed, 4),
                              "MW": round(mw, 2), "LogP": round(logp, 2)})

        scored.sort(key=lambda x: x["QED"], reverse=True)

        return {
            "total_enumerated": len(all_products),
            "lipinski_pass": len(scored),
            "top_k": scored[:top_k]
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    core = "[*:1]c1ccc([*:2])cc1[*:3]"
    rgroups = {1: ["C", "CC", "F"], 2: ["O", "N", "Cl"], 3: ["C", "OC"]}
    result = level_function(core, rgroups, top_k=5)
    if result:
        print(f"Enumerated: {result['total_enumerated']}, Lipinski: {result['lipinski_pass']}")
        for r in result['top_k']:
            print(f"  {r['smiles']}: QED={r['QED']}")
