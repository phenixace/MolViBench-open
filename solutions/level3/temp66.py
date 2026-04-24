from rdkit import Chem


# hERG toxicity risk substructure patterns
HERG_ALERTS = {
    "basic_nitrogen_piperidine": "C1CCNCC1",
    "basic_nitrogen_piperazine": "C1CNCCN1",
    "diphenylmethyl": "C(c1ccccc1)c1ccccc1",
    "long_alkyl_chain": "CCCCCCCC",
    "phenothiazine": "c1ccc2c(c1)Sc1ccccc1N2",
    "biphenyl_amino": "c1ccc(-c2ccccc2)cc1N",
    "basic_amine_aromatic": "c1ccccc1CCN",
    "dimethylamine": "CN(C)",
    "quaternary_nitrogen": "[N+](C)(C)(C)",
    "halogenated_aromatic": "c1cc([F,Cl])ccc1",
}


def level_function(mol):
    """检测分子是否含有 hERG 心脏毒性风险的子结构模式。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        matched_alerts = []
        for name, smarts in HERG_ALERTS.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                continue
            if mol_obj.HasSubstructMatch(pattern):
                matched_alerts.append(name)

        return {
            "num_alerts": len(matched_alerts),
            "alerts": matched_alerts,
            "hERG_risk": len(matched_alerts) > 0
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccc(CCN(C)C)cc1"
    result = level_function(smiles)
    if result:
        print(f"hERG risk: {result['hERG_risk']}")
        print(f"Alerts: {result['alerts']}")
