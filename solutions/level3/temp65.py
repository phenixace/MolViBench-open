from rdkit import Chem


# Brenk structural alerts (subset of commonly used)
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
    "polycyclic_aromatic_gt3": "c1ccc2c(c1)ccc1ccccc12",  # 3+ fused aromatic
    "alpha_halo_ketone": "C(=O)C[F,Cl,Br,I]",
    "mustard": "ClCCN",
}


def level_function(mol):
    """检测分子中是否命中 Brenk 结构警报并列出所有匹配的警报名称。"""
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


if __name__ == "__main__":
    # Test with a molecule containing nitro group
    smiles = "c1ccc([N+](=O)[O-])cc1"
    result = level_function(smiles)
    if result:
        print(f"Brenk alerts: {result['alerts']}")
        print(f"Passes: {result['passes_Brenk']}")
