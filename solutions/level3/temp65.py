from rdkit import Chem

BRENK_ALERTS = {
    "aldehyde": "[CH1](=O)",
    "epoxide": "C1OC1",
    "michael_acceptor": "[#6]=!@[CH2]",
    "peroxide": "OO",
    "sulfonyl_halide": "S(=O)(=O)[F,Cl,Br,I]",
    "acid_halide": "C(=O)[F,Cl,Br,I]",
    "phosphoramide": "NP(=O)(N)N",
    "beta_lactone": "C1(=O)OCC1",
    "aziridine": "C1NC1",
    "isocyanate": "N=C=O",
    "thiocyanate": "SC#N",
    "acyl_cyanide": "C(=O)C#N",
    "azide": "N=[N+]=[N-]",
    "diazo": "[N-]=[N+]=C",
    "disulfide": "SS",
    "thioester": "C(=O)S",
    "hydrazine": "NN",
    "hydroxamic_acid": "ONC=O",
    "nitro": "[N+](=O)[O-]",
    "triflate": "OS(=O)(=O)C(F)(F)F",
    "acetal": "C(O)(O)",
    "enol_ether": "C=CO",
    "alpha_halo_ketone": "C(=O)C[F,Cl,Br,I]",
    "mustard": "ClCCN",
}

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        matched_alerts = []
        for name, smarts in BRENK_ALERTS.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                continue
            if mol_obj.HasSubstructMatch(pattern):
                matched_alerts.append(name)

        return {
            "num_alerts": len(matched_alerts),
            "alerts": matched_alerts,
            "passes_Brenk": len(matched_alerts) == 0
        }
    except Exception as e:
        print(e)
        return None
