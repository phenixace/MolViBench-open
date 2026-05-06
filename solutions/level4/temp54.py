from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import random

def level_function(mol, seed=42):
    try:
        random.seed(seed)
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        current_mol = mol_obj
        rxn_templates = [
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

            random.shuffle(rxn_templates)
            modified = False
            for rxn_sma in rxn_templates:
                rxn = AllChem.ReactionFromSmarts(rxn_sma)
                products = rxn.RunReactants((current_mol,))
                if products:
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
