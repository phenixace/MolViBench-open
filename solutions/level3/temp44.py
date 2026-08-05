from rdkit import Chem
from rdkit.Chem import AllChem


METABOLIC_REACTIONS = [
    ("[CH3:1]>>[CH2:1]O", "Aliphatic oxidation (C-H to C-OH)"),
    ("[cH:1]>>[c:1]O", "Aromatic ring oxidation"),
    ("[NH2:1]>>[NH:1]O", "N-oxidation"),
    ("[SX2:1]>>[S:1](=O)", "S-oxidation"),
    ("[O:1]C>>[O:1]", "O-demethylation"),
    ("[N:1](C)>>[NH:1]", "N-demethylation"),
    ("[C:1](=O)[O:2][C:3]>>[C:1](=O)[OH:2].[OH][C:3]", "Ester hydrolysis"),
    ("[C:1](=O)[C:2]>>[C:1](O)[C:2]", "Ketone reduction"),
]


def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        metabolites = []
        for smarts, name in METABOLIC_REACTIONS:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((mol_obj,))
            for product_set in products:
                for product in product_set:
                    try:
                        Chem.SanitizeMol(product)
                        smi = Chem.MolToSmiles(product)
                        if smi not in [m["smiles"] for m in metabolites]:
                            metabolites.append({"smiles": smi, "reaction": name})
                    except Exception:
                        continue

        return metabolites if metabolites else None
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'CC(=O)Oc1ccccc1OC(C)=O'
    result = level_function(smiles)
    if result:
        for m in result:
            print(f"Output: {m['reaction']}{m['smiles']}")
