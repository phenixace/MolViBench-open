from rdkit import Chem



CYP_SMARTS = {
    "azole": "[nR1]1[cR1][nR1][cR1][cR1]1",
    "thioamide": "[#6]C(=S)N",
    "hydroxamic_acid": "[OH]NC=O",
    "methylenedioxy": "OCOc",
    "furan": "c1ccoc1",
    "thiophene": "c1ccsc1",
    "hydrazine": "NN",
    "isoniazid_like": "c1ccncc1C(=O)NN",
    "terminal_alkyne": "C#[CH]",
    "thiol": "[SH]",
    "primary_amine_aromatic": "c[NH2]",
}


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        matches = {}
        for name, smarts in CYP_SMARTS.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                continue
            has_match = mol_obj.HasSubstructMatch(pattern)
            if has_match:
                matches[name] = True

        is_inhibitor = len(matches) > 0

        return {
            "matched_patterns": list(matches.keys()),
            "num_alerts": len(matches),
            "predicted_CYP_inhibitor": is_inhibitor
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1cnc[nH]1'
    result = level_function(smiles)
    if result:
        print(f"Output: {result['predicted_CYP_inhibitor']}")
        print(f"Output: {result['matched_patterns']}")
