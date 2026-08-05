from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdFMCS
import numpy as np
from collections import defaultdict


def level_function(smiles_list, activities):

    try:
        mols = [(s, Chem.MolFromSmiles(s)) for s in smiles_list]
        mols = [(s, m) for s, m in mols if m is not None]
        if len(mols) < 2:
            return None

        acts = {s: a for (s, _), a in zip(mols, activities)}


        pairs = []
        for i in range(len(mols)):
            for j in range(i + 1, len(mols)):
                smi_i, mol_i = mols[i]
                smi_j, mol_j = mols[j]


                mcs = rdFMCS.FindMCS(
                    [mol_i, mol_j],
                    threshold=0.7,
                    ringMatchesRingOnly=True,
                    completeRingsOnly=True,
                    timeout=5
                )
                if mcs.canceled or mcs.numAtoms < 3:
                    continue


                n_i = mol_i.GetNumHeavyAtoms()
                n_j = mol_j.GetNumHeavyAtoms()
                frac_i = mcs.numAtoms / n_i
                frac_j = mcs.numAtoms / n_j


                if frac_i >= 0.5 and frac_j >= 0.5:
                    core_smarts = mcs.smartsString
                    core_mol = Chem.MolFromSmarts(core_smarts)
                    if core_mol is None:
                        continue


                    match_i = mol_i.GetSubstructMatch(core_mol)
                    match_j = mol_j.GetSubstructMatch(core_mol)
                    if not match_i or not match_j:
                        continue

                    diff_atoms_i = set(range(n_i)) - set(match_i)
                    diff_atoms_j = set(range(n_j)) - set(match_j)

                    act_diff = acts[smi_j] - acts[smi_i]

                    pairs.append({
                        "mol_A": smi_i,
                        "mol_B": smi_j,
                        "core_smarts": core_smarts,
                        "core_atoms": mcs.numAtoms,
                        "diff_atoms_A": len(diff_atoms_i),
                        "diff_atoms_B": len(diff_atoms_j),
                        "activity_A": round(acts[smi_i], 4),
                        "activity_B": round(acts[smi_j], 4),
                        "activity_delta": round(act_diff, 4),
                        "abs_delta": round(abs(act_diff), 4)
                    })


        pairs.sort(key=lambda x: x["abs_delta"], reverse=True)

        return {
            "n_molecules": len(mols),
            "n_mmp_pairs": len(pairs),
            "top_transformations": pairs[:20],
            "max_activity_change": round(pairs[0]["abs_delta"], 4) if pairs else 0.0,
            "mean_activity_change": round(float(np.mean([p["abs_delta"] for p in pairs])), 4) if pairs else 0.0,
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = ['c1ccccc1', 'c1ccc(O)cc1', 'c1ccc(N)cc1', 'c1ccc(F)cc1', 'c1ccc(Cl)cc1', 'c1ccc(OC)cc1']
    acts = [1.0, 2.5, 3.0, 1.2, 1.8, 2.8]
    result = level_function(smiles, acts)
    if result:
        print(f"Output: {result['n_mmp_pairs']}")
        for t in result['top_transformations'][:3]:
            print(f"Output: {t['mol_A']}{t['mol_B']}{t['activity_delta']}")
