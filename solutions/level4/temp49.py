from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None


        pattern = Chem.MolFromSmarts('n')
        has_aromatic_n = mol_obj.HasSubstructMatch(pattern)

        if not has_aromatic_n:
            return None


        rxn = AllChem.ReactionFromSmarts('[n:1]>>[n+:1][O-]')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)


        tpsa = rdMolDescriptors.CalcTPSA(product)

        return {
            "has_aromatic_n": has_aromatic_n,
            "product": product_smiles,
            "tpsa": round(tpsa, 4)
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccncc1'
    print(f'Output: {level_function(smiles)}')
