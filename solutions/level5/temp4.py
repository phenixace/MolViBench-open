from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(pharmacophore_smarts):
    try:
        pattern = Chem.MolFromSmarts(pharmacophore_smarts)
        if pattern is None:
            return None

        fragment_library = [
        ]

        substituent_rxns = [
            '[cH1:1]>>[c:1]O',
            '[cH1:1]>>[c:1]N',
            '[cH1:1]>>[c:1]F',
            '[cH1:1]>>[c:1]Cl',
            '[cH1:1]>>[c:1]C',
            '[cH1:1]>>[c:1]OC',
            '[cH1:1]>>[c:1]C(=O)O',
        ]

        matching_molecules = set()

        for smi in fragment_library:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None and mol.HasSubstructMatch(pattern):
                matching_molecules.add(Chem.MolToSmiles(mol))

        for smi in fragment_library:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            for rxn_sma in substituent_rxns:
                try:
                    rxn = AllChem.ReactionFromSmarts(rxn_sma)
                    products = rxn.RunReactants((mol,))
                    for prod_set in products:
                        for prod in prod_set:
                            try:
                                Chem.SanitizeMol(prod)
                                if prod.HasSubstructMatch(pattern):
                                    prod_smi = Chem.MolToSmiles(prod)
                                    matching_molecules.add(prod_smi)
                            except Exception:
                                pass
                except Exception:
                    pass

        return sorted(list(matching_molecules))
    except Exception as e:
        print(e)
        return None
