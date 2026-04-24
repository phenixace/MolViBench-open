from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


def level_function(pharmacophore_smarts):
    """给定一个已知药效团，生成符合该药效团要求的分子。"""
    try:
        pattern = Chem.MolFromSmarts(pharmacophore_smarts)
        if pattern is None:
            return None

        # Define a library of candidate fragment SMILES to build molecules
        fragment_library = [
            "c1ccccc1",                    # benzene
            "c1ccc(O)cc1",                 # phenol
            "c1ccc(N)cc1",                 # aniline
            "c1ccc(C(=O)O)cc1",            # benzoic acid
            "c1ccc(C(=O)N)cc1",            # benzamide
            "c1ccncc1",                    # pyridine
            "c1ccc2[nH]ccc2c1",            # indole
            "c1ccc(NC(=O)C)cc1",           # acetanilide
            "CC(=O)Oc1ccccc1C(=O)O",       # aspirin
            "CC(C)Cc1ccc(C(C)C(=O)O)cc1",  # ibuprofen
            "c1ccc(CC(=O)O)cc1",           # phenylacetic acid
            "c1ccc(CCN)cc1",               # phenethylamine
            "c1ccc(CCNC)cc1",              # N-methyl phenethylamine
            "c1ccc(CCC(=O)O)cc1",          # hydrocinnamic acid
            "OC(=O)c1cccc(O)c1",           # 3-hydroxybenzoic acid
            "Nc1ccc(O)cc1",                # aminophenol
            "c1ccc(NS(=O)(=O)c2ccccc2)cc1",  # sulfonamide
            "c1cnc2ccccc2c1",              # quinoline
            "c1ccc2c(c1)cccc2O",           # naphthol
            "CC(=O)Nc1ccc(O)cc1",          # paracetamol
            "c1ccc(C(O)c2ccccc2)cc1",      # benzhydrol
            "c1ccc(Oc2ccccc2)cc1",         # diphenyl ether
            "c1ccc(-c2ccccn2)cc1",         # 2-phenylpyridine
            "c1cc(O)c(O)cc1CC=C",          # eugenol-like
            "OC(=O)CCc1ccc(O)c(O)c1",      # caffeic acid analog
        ]

        # Additionally modify fragments
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

        # Check base fragments
        for smi in fragment_library:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None and mol.HasSubstructMatch(pattern):
                matching_molecules.add(Chem.MolToSmiles(mol))

        # Generate modified fragments and check
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


if __name__ == "__main__":
    # Search for molecules containing a phenol substructure
    result = level_function("[OH]c1ccccc1")
    print(f"result: {result}")
