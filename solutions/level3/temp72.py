from rdkit import Chem

def level_function(generated_smiles, reference_smiles=None):
    try:
        if not generated_smiles:
            return None

        total = len(generated_smiles)

        valid_mols = []
        valid_smiles_set = set()
        for smi in generated_smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                canonical = Chem.MolToSmiles(mol)
                valid_mols.append(canonical)
                valid_smiles_set.add(canonical)

        validity = len(valid_mols) / total if total > 0 else 0.0

        uniqueness = len(valid_smiles_set) / len(valid_mols) if valid_mols else 0.0

        novelty = 1.0
        if reference_smiles:
            ref_set = set()
            for smi in reference_smiles:
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    ref_set.add(Chem.MolToSmiles(mol))

            novel_count = sum(1 for s in valid_smiles_set if s not in ref_set)
            novelty = novel_count / len(valid_smiles_set) if valid_smiles_set else 0.0

        return {
            "total": total,
            "valid": len(valid_mols),
            "unique": len(valid_smiles_set),
            "validity": round(validity, 4),
            "uniqueness": round(uniqueness, 4),
            "novelty": round(novelty, 4)
        }
    except Exception as e:
        print(e)
        return None
