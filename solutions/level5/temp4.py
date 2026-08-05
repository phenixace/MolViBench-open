from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


def level_function(pharmacophore_smarts):

    try:
        pattern = Chem.MolFromSmarts(pharmacophore_smarts)
        if pattern is None:
            return None


        fragment_library = [
            "c1ccccc1",
            "c1ccc(O)cc1",
            "c1ccc(N)cc1",
            "c1ccc(C(=O)O)cc1",
            "c1ccc(C(=O)N)cc1",
            "c1ccncc1",
            "c1ccc2[nH]ccc2c1",
            "c1ccc(NC(=O)C)cc1",
            "CC(=O)Oc1ccccc1C(=O)O",
            "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
            "c1ccc(CC(=O)O)cc1",
            "c1ccc(CCN)cc1",
            "c1ccc(CCNC)cc1",
            "c1ccc(CCC(=O)O)cc1",
            "OC(=O)c1cccc(O)c1",
            "Nc1ccc(O)cc1",
            "c1ccc(NS(=O)(=O)c2ccccc2)cc1",
            "c1cnc2ccccc2c1",
            "c1ccc2c(c1)cccc2O",
            "CC(=O)Nc1ccc(O)cc1",
            "c1ccc(C(O)c2ccccc2)cc1",
            "c1ccc(Oc2ccccc2)cc1",
            "c1ccc(-c2ccccn2)cc1",
            "c1cc(O)c(O)cc1CC=C",
            "OC(=O)CCc1ccc(O)c(O)c1",
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


if __name__ == '__main__':
    result = level_function('[OH]c1ccccc1')
    print(f'Output: {result}')
