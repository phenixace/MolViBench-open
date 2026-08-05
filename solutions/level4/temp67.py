from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        mw = Descriptors.MolWt(mol_obj)
        logp = Crippen.MolLogP(mol_obj)
        hbd = rdMolDescriptors.CalcNumHBD(mol_obj)
        hba = rdMolDescriptors.CalcNumHBA(mol_obj)

        rules = [
            {"rule": "MW < 500", "value": round(mw, 2), "passes": mw < 500},
            {"rule": "LogP < 5", "value": round(logp, 4), "passes": logp < 5},
            {"rule": "HBD ≤ 5", "value": hbd, "passes": hbd <= 5},
            {"rule": "HBA ≤ 10", "value": hba, "passes": hba <= 10},
        ]

        violations = sum(1 for r in rules if not r["passes"])

        return {
            "smiles": Chem.MolToSmiles(mol_obj),
            "rules": rules,
            "violations": violations,
            "passes_Lipinski": violations <= 1
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'CC(C)Cc1ccc(C(C)C(=O)O)cc1'
    result = level_function(smiles)
    if result:
        for r in result['rules']:
            status = '✓' if r['passes'] else '✗'
            print(f"Output: {status}{r['rule']}{r['value']}")
        print(f"Output: {result['passes_Lipinski']}")
