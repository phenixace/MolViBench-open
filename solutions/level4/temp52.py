from rdkit import Chem
from rdkit.Chem import AllChem, Crippen


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        logp = Crippen.MolLogP(mol_obj)
        original_logp = logp
        action = "none"
        final_mol = mol_obj

        if logp > 3:

            rxn = AllChem.ReactionFromSmarts('[cH:1]>>[c:1]O')
            products = rxn.RunReactants((mol_obj,))
            if products:
                final_mol = products[0][0]
                Chem.SanitizeMol(final_mol)
                action = "add_OH"
            else:
                rxn2 = AllChem.ReactionFromSmarts('[CH3:1]>>[CH2:1]O')
                products = rxn2.RunReactants((mol_obj,))
                if products:
                    final_mol = products[0][0]
                    Chem.SanitizeMol(final_mol)
                    action = "add_OH"
        elif logp < 0:

            rxn = AllChem.ReactionFromSmarts('[cH:1]>>[c:1]C')
            products = rxn.RunReactants((mol_obj,))
            if products:
                final_mol = products[0][0]
                Chem.SanitizeMol(final_mol)
                action = "add_CH3"
            else:
                rxn2 = AllChem.ReactionFromSmarts('[NH2:1]>>[NH:1]C')
                products = rxn2.RunReactants((mol_obj,))
                if products:
                    final_mol = products[0][0]
                    Chem.SanitizeMol(final_mol)
                    action = "add_CH3"

        new_logp = Crippen.MolLogP(final_mol)
        in_range = 1.0 <= new_logp <= 3.0

        return {
            "original_LogP": round(original_logp, 4),
            "action": action,
            "final_smiles": Chem.MolToSmiles(final_mol),
            "new_LogP": round(new_logp, 4),
            "in_range_1_3": in_range
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccc(CCCC)cc1'
    print(f'Output: {level_function(smiles)}')
