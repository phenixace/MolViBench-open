from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


def level_function(target_mw_min=200, target_mw_max=250):

    try:
        start_mol = Chem.MolFromSmiles("c1ccccc1")

        substituent_rxns = [
            ('[cH:1]>>[c:1]C', 'methyl'),
            ('[cH:1]>>[c:1]O', 'OH'),
            ('[cH:1]>>[c:1]N', 'NH2'),
            ('[cH:1]>>[c:1]F', 'F'),
            ('[cH:1]>>[c:1]Cl', 'Cl'),
            ('[cH:1]>>[c:1]CC', 'ethyl'),
            ('[cH:1]>>[c:1]OC', 'OMe'),
            ('[cH:1]>>[c:1]C(=O)O', 'COOH'),
            ('[cH:1]>>[c:1]C(=O)N', 'CONH2'),
            ('[cH:1]>>[c:1]Br', 'Br'),
        ]

        def backtrack(mol, depth, path, max_depth=8):
            mw = Descriptors.MolWt(mol)
            if target_mw_min <= mw <= target_mw_max:
                return Chem.MolToSmiles(mol), path, mw
            if mw > target_mw_max or depth >= max_depth:
                return None, None, None

            for rxn_sma, name in substituent_rxns:
                rxn = AllChem.ReactionFromSmarts(rxn_sma)
                products = rxn.RunReactants((mol,))
                if products:
                    prod = products[0][0]
                    try:
                        Chem.SanitizeMol(prod)
                        result, result_path, result_mw = backtrack(
                            prod, depth + 1, path + [name])
                        if result is not None:
                            return result, result_path, result_mw
                    except Exception:
                        continue
            return None, None, None

        result_smi, result_path, result_mw = backtrack(start_mol, 0, [])

        if result_smi:
            return {
                "start": "c1ccccc1",
                "final_smiles": result_smi,
                "final_MW": round(result_mw, 2),
                "steps": result_path,
                "num_steps": len(result_path)
            }
        return None
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    result = level_function(200, 250)
    if result:
        print(f"Output: {result['final_smiles']}{result['final_MW']}")
        print(f"Output: {result['steps']}")
