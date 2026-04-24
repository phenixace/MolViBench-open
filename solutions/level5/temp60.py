import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Crippen, BRICS, FilterCatalog
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker


def level_function(scaffold_smiles, rgroup_lists, n_select=50, seed=42):
    """组合库设计完整流程：给定骨架和 R-group 列表 → 枚举所有组合 → 级联过滤（Lipinski + Veber + PAINS + Brenk）→ MaxMin 多样性选择 Top-50 → 输出候选及性质摘要。"""
    try:
        np.random.seed(seed)

        scaffold = Chem.MolFromSmiles(scaffold_smiles)
        if scaffold is None:
            return None

        # Enumerate combinations using ReplaceSubstructs
        # rgroup_lists: list of lists, each inner list = SMILES for one R-group position
        # Use [*:1], [*:2], etc. as attachment points

        # Simple enumeration: for each R-group position, replace dummy atoms
        products = set()
        if len(rgroup_lists) == 0:
            return None

        # Start with scaffold, iteratively replace each R-group
        current_mols = [scaffold_smiles]
        for pos_idx, rgroups in enumerate(rgroup_lists):
            next_mols = []
            dummy_smarts = f"[#{0}]"  # Use wildcard atom [*]
            for base_smi in current_mols:
                base_mol = Chem.MolFromSmiles(base_smi)
                if base_mol is None:
                    continue
                for rg_smi in rgroups:
                    rg_mol = Chem.MolFromSmiles(rg_smi)
                    if rg_mol is None:
                        continue
                    try:
                        # Direct combination via SMILES concatenation at attachment point
                        combined_smi = base_smi.replace(f"[*:{pos_idx+1}]", rg_smi, 1)
                        combined_mol = Chem.MolFromSmiles(combined_smi)
                        if combined_mol is not None:
                            Chem.SanitizeMol(combined_mol)
                            next_mols.append(Chem.MolToSmiles(combined_mol))
                    except Exception:
                        continue
            if next_mols:
                current_mols = list(set(next_mols))

        products = set(current_mols)
        # Remove the original scaffold if it's still there
        products.discard(scaffold_smiles)

        n_enumerated = len(products)

        # Parse all products
        parsed = []
        for smi in products:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                parsed.append((smi, mol))

        # Cascade filtering
        # 1. Lipinski
        lipinski_pass = []
        for smi, mol in parsed:
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            if mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10:
                lipinski_pass.append((smi, mol))
        n_lipinski = len(lipinski_pass)

        # 2. Veber
        veber_pass = []
        for smi, mol in lipinski_pass:
            tpsa = Descriptors.TPSA(mol)
            rotb = Descriptors.NumRotatableBonds(mol)
            if tpsa <= 140 and rotb <= 10:
                veber_pass.append((smi, mol))
        n_veber = len(veber_pass)

        # 3. PAINS filter
        pains_params = FilterCatalog.FilterCatalogParams()
        pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        pains_catalog = FilterCatalog.FilterCatalog(pains_params)

        pains_pass = []
        for smi, mol in veber_pass:
            if pains_catalog.GetFirstMatch(mol) is None:
                pains_pass.append((smi, mol))
        n_pains = len(pains_pass)

        # 4. Brenk filter
        brenk_params = FilterCatalog.FilterCatalogParams()
        brenk_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
        brenk_catalog = FilterCatalog.FilterCatalog(brenk_params)

        brenk_pass = []
        for smi, mol in pains_pass:
            if brenk_catalog.GetFirstMatch(mol) is None:
                brenk_pass.append((smi, mol))
        n_brenk = len(brenk_pass)

        # MaxMin diversity selection
        if len(brenk_pass) <= n_select:
            selected = brenk_pass
        else:
            fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for _, m in brenk_pass]

            def dist_fn(i, j):
                return 1.0 - DataStructs.TanimotoSimilarity(fps[i], fps[j])

            picker = MaxMinPicker()
            picks = picker.LazyPick(dist_fn, len(fps), n_select, seed=seed)
            selected = [brenk_pass[i] for i in picks]

        # Output properties summary
        results = []
        for smi, mol in selected:
            results.append({
                "smiles": smi,
                "MW": round(Descriptors.MolWt(mol), 2),
                "LogP": round(Crippen.MolLogP(mol), 4),
                "TPSA": round(Descriptors.TPSA(mol), 2),
                "QED": round(Descriptors.qed(mol), 4),
                "HBD": Descriptors.NumHDonors(mol),
                "HBA": Descriptors.NumHAcceptors(mol),
            })

        return {
            "scaffold": scaffold_smiles,
            "n_rgroup_positions": len(rgroup_lists),
            "pipeline_summary": {
                "enumerated": n_enumerated,
                "after_lipinski": n_lipinski,
                "after_veber": n_veber,
                "after_PAINS": n_pains,
                "after_Brenk": n_brenk,
                "final_selected": len(results),
            },
            "selected_molecules": results,
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    scaffold = "c1ccc([*:1])c([*:2])c1"
    rgroups = [
        ["F", "Cl", "Br", "O", "N", "C"],
        ["C(=O)O", "C(=O)N", "C#N", "C(F)(F)F", "OC"],
    ]
    result = level_function(scaffold, rgroups, n_select=10)
    if result:
        print(f"Pipeline: {result['pipeline_summary']}")
        for m in result['selected_molecules'][:5]:
            print(f"  {m['smiles']} (QED={m['QED']})")
