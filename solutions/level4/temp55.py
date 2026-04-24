from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import random


def level_function(mol, seed=42):
    """给定分子 → 迭代替换不同侧链 → 每轮计算 QED → 当 QED 连续两轮不再提升时停止 → 输出最优分子及其 QED。"""
    try:
        random.seed(seed)
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        substituents = [
            ('[cH:1]>>[c:1]C', 'add_methyl'),
            ('[cH:1]>>[c:1]O', 'add_OH'),
            ('[cH:1]>>[c:1]N', 'add_NH2'),
            ('[cH:1]>>[c:1]F', 'add_F'),
            ('[cH:1]>>[c:1]Cl', 'add_Cl'),
            ('[cH:1]>>[c:1]OC', 'add_OMe'),
            ('[cH:1]>>[c:1]C(=O)N', 'add_amide'),
        ]

        current_mol = mol_obj
        best_qed = Descriptors.qed(current_mol)
        best_mol = current_mol
        prev_qed = -1
        no_improve_count = 0
        max_iters = 20
        history = []

        for iteration in range(max_iters):
            current_qed = Descriptors.qed(current_mol)
            history.append({
                "iter": iteration,
                "QED": round(current_qed, 4),
                "smiles": Chem.MolToSmiles(current_mol)
            })

            if current_qed > best_qed:
                best_qed = current_qed
                best_mol = current_mol
                no_improve_count = 0
            else:
                no_improve_count += 1

            if no_improve_count >= 2:
                break

            # Try all substituents, pick best
            best_next = None
            best_next_qed = current_qed
            random.shuffle(substituents)
            for rxn_sma, name in substituents:
                rxn = AllChem.ReactionFromSmarts(rxn_sma)
                products = rxn.RunReactants((current_mol,))
                for prod_set in products:
                    try:
                        prod = prod_set[0]
                        Chem.SanitizeMol(prod)
                        q = Descriptors.qed(prod)
                        if q > best_next_qed:
                            best_next_qed = q
                            best_next = prod
                    except Exception:
                        continue

            if best_next is not None:
                current_mol = best_next
            else:
                break

            prev_qed = current_qed

        return {
            "best_smiles": Chem.MolToSmiles(best_mol),
            "best_QED": round(best_qed, 4),
            "iterations": len(history),
            "history": history
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccccc1"
    result = level_function(smiles)
    if result:
        print(f"Best: {result['best_smiles']}, QED: {result['best_QED']}")
