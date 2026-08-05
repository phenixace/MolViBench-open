from rdkit import Chem


TOXIC_SUBSTRUCTURES = {
    "Aromatic nitro compound": "[$(a[N+](=O)[O-]),$(a[N](=O)=O)]",
    "Aromatic amine": "[NH2]a",
    "Aldehyde": "[CH1](=O)",
    "Epoxide": "C1OC1",
    "Acyl halide": "[CX3](=[OX1])[F,Cl,Br,I]",
    "Isocyanate": "[N]=[C]=[O]",
    "Azo compound": "[N]=[N]",
    "Nitroso group": "[N]=O",
    "Sulfonate ester": "S(=O)(=O)O[C,c]",
    "Phosphate ester": "P(=O)(O)(O)O",
    "Peroxide": "OO",
    "Michael acceptor": "[CH2]=[CH][C,S,N](=O)",
    "Polyhalogenated hydrocarbon": "[CX4]([F,Cl,Br,I])([F,Cl,Br,I])",
}


def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        found_alerts = []
        for name, smarts in TOXIC_SUBSTRUCTURES.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                continue
            if mol_obj.HasSubstructMatch(pattern):
                found_alerts.append(name)

        return found_alerts if found_alerts else []
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccc([N+](=O)[O-])cc1N'
    print(f'Output: {level_function(smiles)}')
