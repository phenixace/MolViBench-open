from rdkit import Chem

def level_function(smiles_list):
    try:
        results = []
        success_count = 0
        fail_count = 0
        fixed_count = 0

        for smi in smiles_list:
            original = smi
            mol = Chem.MolFromSmiles(smi)

            if mol is not None:
                success_count += 1
                results.append({
                    "original": original,
                    "status": "valid",
                    "canonical": Chem.MolToSmiles(mol)
                })
                continue

            mol = Chem.MolFromSmiles(fixed_smi)
            if mol is not None:
                fixed_count += 1
                success_count += 1
                results.append({
                    "original": original,
                    "status": "fixed_whitespace",
                    "canonical": Chem.MolToSmiles(mol)
                })
                continue

            open_count = fixed_smi.count('(')
            close_count = fixed_smi.count(')')
            if open_count > close_count:
                fixed_smi += ')' * (open_count - close_count)
                mol = Chem.MolFromSmiles(fixed_smi)
                if mol is not None:
                    fixed_count += 1
                    success_count += 1
                    results.append({
                        "original": original,
                        "status": "fixed_parentheses",
                        "canonical": Chem.MolToSmiles(mol)
                    })
                    continue

            if close_count > open_count:
                fixed_smi2 = fixed_smi[::-1].replace(')', '', close_count - open_count)[::-1]
                mol = Chem.MolFromSmiles(fixed_smi2)
                if mol is not None:
                    fixed_count += 1
                    success_count += 1
                    results.append({
                        "original": original,
                        "status": "fixed_parentheses",
                        "canonical": Chem.MolToSmiles(mol)
                    })
                    continue

            fail_count += 1
            results.append({
                "original": original,
                "status": "failed",
                "canonical": None
            })

        return {
            "total": len(smiles_list),
            "success": success_count,
            "fixed": fixed_count,
            "failed": fail_count,
            "details": results
        }
    except Exception as e:
        print(e)
        return None
