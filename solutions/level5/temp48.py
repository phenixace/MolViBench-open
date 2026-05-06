from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

PHARMACOPHORE_SMARTS = {
    "H-bond donor": "[#7H,#8H]",
    "H-bond acceptor": "[#7,#8]",
    "Aromatic ring": "a1aaaaa1",
    "Hydrophobic": "[CH2,CH3]",
}

def level_function(pharmacophore_type="Aromatic ring"):
    try:
        if pharmacophore_type not in PHARMACOPHORE_SMARTS:
            return None

        smarts = PHARMACOPHORE_SMARTS[pharmacophore_type]

        base_molecules = [
            "c1ccccc1", "c1ccncc1", "c1ccoc1", "c1ccsc1",
            "c1ccc(O)cc1", "c1ccc(N)cc1", "c1ccc(F)cc1",
            "c1ccc(C)cc1", "c1ccc(CC)cc1", "c1ccc2ccccc2c1",
        ]

        rxns = [
            '[cH:1]>>[c:1]C',
            '[cH:1]>>[c:1]O',
            '[cH:1]>>[c:1]F',
            '[cH:1]>>[c:1]N',
        ]

        candidates = set()
        for base_smi in base_molecules:
            base = Chem.MolFromSmiles(base_smi)
            if base is None:
                continue
            pattern = Chem.MolFromSmarts(smarts)
            if not base.HasSubstructMatch(pattern):
                continue
            candidates.add(base_smi)
            for rxn_smarts in rxns:
                rxn = AllChem.ReactionFromSmarts(rxn_smarts)
                products = rxn.RunReactants((base,))
                for ps in products:
                    for p in ps:
                        try:
                            Chem.SanitizeMol(p)
                            smi = Chem.MolToSmiles(p)
                            candidates.add(smi)
                        except Exception:
                            continue

        filtered = []
        for smi in candidates:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            if tpsa < 90:
                mw = Descriptors.MolWt(mol)
                filtered.append({
                    'smiles': smi,
                    'tpsa': round(tpsa, 2),
                    'mw': round(mw, 2)
                })

        filtered.sort(key=lambda x: x['mw'])
        return filtered[:3]
    except Exception as e:
        print(e)
        return None
