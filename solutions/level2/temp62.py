from rdkit import Chem
from rdkit.Chem import AllChem, rdChemReactions


def level_function(mol, reaction_smarts):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        if rxn is None:
            return None


        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None


        result = []
        for prod_set in products[:1]:
            for prod in prod_set:
                Chem.SanitizeMol(prod)
                mapping = {}
                for atom in prod.GetAtoms():
                    map_num = atom.GetAtomMapNum()
                    if map_num > 0:
                        mapping[atom.GetIdx()] = map_num
                result.append({
                    "product_smiles": Chem.MolToSmiles(prod),
                    "atom_mapping": mapping
                })
        return result
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccccc1O'
    rxn_smarts = '[OH:1]>>[Cl:1]'
    print(f'Output: {level_function(smiles, rxn_smarts)}')
