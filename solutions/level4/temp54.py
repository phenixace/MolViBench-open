from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import random


def level_function(mol, seed=42):
    """给定分子 → 反复在随机位置添加极性基团 → 每轮计算 TPSA → 直到 TPSA 落入 [60, 90] 区间 → 输出最终分子和迭代次数。"""
    try:
        random.seed(seed)
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        current_mol = mol_obj
        polar_groups = ['O', 'N', 'O', 'N']  # OH, NH2 additions
        rxn_templates = [
            '[cH:1]>>[c:1]O',     # Add -OH to aromatic
            '[cH:1]>>[c:1]N',     # Add -NH2 to aromatic
            '[CH3:1]>>[CH2:1]O',  # Add -OH to methyl
            '[CH2:1]>>[CH:1]O',   # Add -OH to methylene
        ]

        iterations = 0
        max_iterations = 20
        history = []

        while iterations < max_iterations:
            tpsa = Descriptors.TPSA(current_mol)
            history.append({"iter": iterations, "TPSA": round(tpsa, 2),
                           "smiles": Chem.MolToSmiles(current_mol)})

            if 60 <= tpsa <= 90:
                return {
                    "final_smiles": Chem.MolToSmiles(current_mol),
                    "final_TPSA": round(tpsa, 2),
                    "iterations": iterations,
                    "history": history
                }

            # Try adding a polar group
            random.shuffle(rxn_templates)
            modified = False
            for rxn_sma in rxn_templates:
                rxn = AllChem.ReactionFromSmarts(rxn_sma)
                products = rxn.RunReactants((current_mol,))
                if products:
                    # Pick a random product
                    prod_idx = random.randint(0, len(products) - 1)
                    prod = products[prod_idx][0]
                    try:
                        Chem.SanitizeMol(prod)
                        current_mol = prod
                        modified = True
                        break
                    except Exception:
                        continue

            if not modified:
                break
            iterations += 1

        tpsa = Descriptors.TPSA(current_mol)
        return {
            "final_smiles": Chem.MolToSmiles(current_mol),
            "final_TPSA": round(tpsa, 2),
            "iterations": iterations,
            "converged": 60 <= tpsa <= 90,
            "history": history
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccccc1"  # Benzene, TPSA = 0
    result = level_function(smiles)
    if result:
        print(f"Final: {result['final_smiles']}, TPSA: {result['final_TPSA']}, Iters: {result['iterations']}")
