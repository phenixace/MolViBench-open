from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen
from rdkit.Chem.Scaffolds import MurckoScaffold


def level_function(mol_smiles, replacement_rings=None):

    try:
        mol = Chem.MolFromSmiles(mol_smiles)
        if mol is None:
            return None


        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold_smi = Chem.MolToSmiles(scaffold)


        if replacement_rings is None:
            replacement_rings = [
                "c1ccncc1",
                "c1ccoc1",
                "c1ccsc1",
                "c1cc[nH]c1",
                "c1cnc2ccccc2n1",
                "C1CCNCC1",
                "C1CCOCC1",
                "c1ccc2[nH]ccc2c1",
                "c1cnc[nH]1",
                "c1ccnnc1",
            ]


        orig_logp = Crippen.MolLogP(mol)
        orig_tpsa = Descriptors.TPSA(mol)
        orig_qed = Descriptors.qed(mol)

        results = []
        for ring_smi in replacement_rings:
            ring_mol = Chem.MolFromSmiles(ring_smi)
            if ring_mol is None:
                continue

            replaced = AllChem.ReplaceSubstructs(mol, scaffold, ring_mol)
            if not replaced:
                continue

            for prod in replaced:
                try:
                    Chem.SanitizeMol(prod)
                    new_smi = Chem.MolToSmiles(prod)
                    new_logp = Crippen.MolLogP(prod)
                    new_tpsa = Descriptors.TPSA(prod)
                    new_qed = Descriptors.qed(prod)

                    results.append({
                        "new_ring": ring_smi,
                        "morphed_smiles": new_smi,
                        "LogP": round(new_logp, 4),
                        "TPSA": round(new_tpsa, 2),
                        "QED": round(new_qed, 4),
                        "delta_LogP": round(new_logp - orig_logp, 4),
                        "delta_TPSA": round(new_tpsa - orig_tpsa, 2),
                        "delta_QED": round(new_qed - orig_qed, 4)
                    })
                    break
                except Exception:
                    continue

        return {
            "original": {
                "smiles": mol_smiles,
                "scaffold": scaffold_smi,
                "LogP": round(orig_logp, 4),
                "TPSA": round(orig_tpsa, 2),
                "QED": round(orig_qed, 4)
            },
            "morphed": results
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    result = level_function('c1ccc(NC(=O)C)cc1')
    if result:
        print(f"Output: {result['original']['scaffold']}")
        for m in result['morphed'][:3]:
            print(f"Output: {m['new_ring']}{m['morphed_smiles']}{m['delta_QED']}")
