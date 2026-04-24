from rdkit import Chem


def level_function(smiles_list):
    """给定一组可能包含错误的 SMILES → 逐个尝试解析 → 对无效 SMILES 尝试修复（去空格、补全括号）→ 统计成功解析和失败的数量。"""
    try:
        results = []
        success_count = 0
        fail_count = 0
        fixed_count = 0

        for smi in smiles_list:
            original = smi
            # First try direct parse
            mol = Chem.MolFromSmiles(smi)

            if mol is not None:
                success_count += 1
                results.append({
                    "original": original,
                    "status": "valid",
                    "canonical": Chem.MolToSmiles(mol)
                })
                continue

            # Attempt fixes
            fixed_smi = smi.strip()  # Remove whitespace
            fixed_smi = fixed_smi.replace(" ", "")  # Remove internal spaces

            # Try after whitespace fix
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

            # Try adding missing closing parentheses
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

            # Try removing closing parentheses
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

            # Failed to fix
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


if __name__ == "__main__":
    smiles = ["CCO", "  c1ccccc1  ", "CC(=O", "INVALID", "CC(O)C", "c1ccc(cc1"]
    result = level_function(smiles)
    if result:
        print(f"Success: {result['success']}, Fixed: {result['fixed']}, Failed: {result['failed']}")
