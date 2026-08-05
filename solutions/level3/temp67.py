from rdkit import Chem



AMES_ALERTS = {
    "aromatic_amine": "[NH2]c",
    "nitroso": "N=O",
    "nitro_aromatic": "c[N+](=O)[O-]",
    "aromatic_azo": "c/N=N/c",
    "alkyl_halide": "[CH2][Cl,Br,I]",
    "epoxide": "C1OC1",
    "aromatic_nitro": "[cR1][N+](=O)[O-]",
    "polycyclic_aromatic": "c1ccc2c(c1)ccc1ccccc12",
    "hydrazine": "[NH]N",
    "azide": "[N-]=[N+]=[N-]",
    "hydroxylamine": "[NH]O",
    "mustard": "N(CCCl)CCCl",
    "aminoazo": "Nc1ccc(N=Nc2ccccc2)cc1",
    "Michael_acceptor_ketone": "C(=O)C=C",
    "acyl_hydrazide": "C(=O)NN",
}


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        matched_alerts = []
        for name, smarts in AMES_ALERTS.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                continue
            if mol_obj.HasSubstructMatch(pattern):
                matched_alerts.append(name)

        return {
            "num_alerts": len(matched_alerts),
            "alerts": matched_alerts,
            "Ames_risk": len(matched_alerts) > 0
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccc(N)cc1'
    result = level_function(smiles)
    if result:
        print(f"Output: {result['Ames_risk']}")
        print(f"Output: {result['alerts']}")
