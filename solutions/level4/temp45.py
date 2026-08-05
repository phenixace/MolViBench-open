from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None


        pattern = Chem.MolFromSmarts('[C:1](=O)[O:2][C:3]')
        has_ester = mol_obj.HasSubstructMatch(pattern)

        if not has_ester:
            return None


        rxn = AllChem.ReactionFromSmarts(
            '[C:1](=O)[O:2][C:3]>>[C:1](=O)O.[C:3]O'
        )
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None


        fragments = []
        for product in products[0]:
            Chem.SanitizeMol(product)
            smi = Chem.MolToSmiles(product)
            mw = rdMolDescriptors.CalcExactMolWt(product)
            fragments.append({
                "smiles": smi,
                "molecular_weight": round(mw, 4)
            })

        return {
            "has_ester": has_ester,
            "fragments": fragments
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'CC(=O)OCC'
    print(f'Output: {level_function(smiles)}')
