from rdkit import Chem



BBB_RULES = {
    "MW_lt_400": lambda mol: Chem.Descriptors.MolWt(mol) < 400,
    "LogP_1_to_3": lambda mol: 1.0 <= Chem.Crippen.MolLogP(mol) <= 3.0,
    "HBD_lt_3": lambda mol: Chem.rdMolDescriptors.CalcNumHBD(mol) < 3,
    "HBA_lt_7": lambda mol: Chem.rdMolDescriptors.CalcNumHBA(mol) < 7,
    "TPSA_lt_90": lambda mol: Chem.Descriptors.TPSA(mol) < 90,
    "RotBonds_lt_8": lambda mol: Chem.rdMolDescriptors.CalcNumRotatableBonds(mol) < 8,
}


def level_function(mol):

    try:
        from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors

        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        mw = Descriptors.MolWt(mol_obj)
        logp = Crippen.MolLogP(mol_obj)
        hbd = rdMolDescriptors.CalcNumHBD(mol_obj)
        hba = rdMolDescriptors.CalcNumHBA(mol_obj)
        tpsa = Descriptors.TPSA(mol_obj)
        rotbonds = rdMolDescriptors.CalcNumRotatableBonds(mol_obj)

        rules = {
            "MW < 400": mw < 400,
            "1 ≤ LogP ≤ 3": 1.0 <= logp <= 3.0,
            "HBD < 3": hbd < 3,
            "HBA < 7": hba < 7,
            "TPSA < 90": tpsa < 90,
            "RotBonds < 8": rotbonds < 8,
        }

        passes_all = all(rules.values())
        num_passed = sum(rules.values())

        return {
            "properties": {
                "MW": round(mw, 2), "LogP": round(logp, 2),
                "HBD": hbd, "HBA": hba,
                "TPSA": round(tpsa, 2), "RotBonds": rotbonds
            },
            "rules": rules,
            "BBB_permeable": passes_all,
            "rules_passed": f"{num_passed}/{len(rules)}"
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccc(NC(=O)C)cc1'
    result = level_function(smiles)
    if result:
        print(f"Output: {result['BBB_permeable']}")
        print(f"Output: {result['rules_passed']}")
