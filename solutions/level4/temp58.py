from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, DataStructs, rdMolDescriptors
import random

def level_function(mol_smi, iterations=15, seed=42):
    random.seed(seed)
    current_mol = Chem.MolFromSmiles(mol_smi)
    if not current_mol: return None

    orig_fp = AllChem.GetMorganFingerprintAsBitVect(current_mol, 2, nBits=2048)

    replacements = [
        ('[CX4H3:1]>>[*:1]O', 'Hydroxyl'),
        ('[CX4H3:1]>>[*:1]C(=O)O', 'Carboxyl'),
        ('[CX4H3:1]>>[*:1]C(=O)N', 'Amide'),
        ('[CX4H3:1]>>[*:1]S(=O)(=O)C', 'Methylsulfonyl'),
        ('[CX4H2:1]>>[*:1][NH1]', 'Introduce an amino group into the chain'),
        ('[CX4H2:1]>>[*:1]O', 'Introduce an ether linkage into the chain'),
        ('[H:1][c:2]>>[c:2]F', 'Fluorination (increase metabolic stability)')
    ]

    def get_detailed_score(m):
        m.UpdatePropertyCache(strict=False)
        mw = Descriptors.MolWt(m)
        logp = Crippen.MolLogP(m)
        violations = sum([mw >= 500, logp >= 5])
        score = (violations * 100) + abs(logp - 2.0)
        return score, {"MW": round(mw, 2), "LogP": round(logp, 2), "Violations": violations}

    current_score, current_props = get_detailed_score(current_mol)

    for i in range(iterations):
        rxn_smarts, label = random.choice(replacements)
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)

        products = rxn.RunReactants((current_mol,))
        if not products: continue

        test_mol = random.choice(products)[0]

        try:
            Chem.SanitizeMol(test_mol)

            new_fp = AllChem.GetMorganFingerprintAsBitVect(test_mol, 2, nBits=2048)
            sim = DataStructs.TanimotoSimilarity(orig_fp, new_fp)

            if sim < 0.5: continue

            new_score, new_props = get_detailed_score(test_mol)

            if new_score < current_score:
                current_mol = test_mol
                current_score = new_score
                current_props = new_props
        except:
            continue

    return {
        "final_smiles": Chem.MolToSmiles(current_mol),
        "final_props": current_props,
        "similarity": round(DataStructs.TanimotoSimilarity(orig_fp, AllChem.GetMorganFingerprintAsBitVect(current_mol, 2)), 4)
    }

if __name__ == '__main__':
    input_smi = 'c1ccc(CCCCCCC)cc1'
    result = level_function(input_smi)
    print(f'Output: {result}')
