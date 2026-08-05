from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None


        num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol_obj)
        has_aromatic_ring = num_aromatic_rings > 0

        if not has_aromatic_ring:
            return None


        rxn = AllChem.ReactionFromSmarts('[cH:1]>>[c:1]Br')

        products1 = rxn.RunReactants((mol_obj,))
        if not products1:
            return None


        all_disubstituted = set()
        for prod_tuple in products1:
            p1 = prod_tuple[0]
            try:
                Chem.SanitizeMol(p1)
            except Exception:
                continue

            products2 = rxn.RunReactants((p1,))
            if products2:
                for prod_tuple2 in products2:
                    p2 = prod_tuple2[0]
                    try:
                        Chem.SanitizeMol(p2)
                        smi = Chem.MolToSmiles(p2)
                        all_disubstituted.add(smi)
                    except Exception:
                        continue

        num_products = len(all_disubstituted)
        product_list = sorted(list(all_disubstituted))

        return {
            "has_aromatic_ring": has_aromatic_ring,
            "products": product_list,
            "num_products": num_products
        }
    except Exception as e:
        print(e)
        return None

if __name__ == '__main__':
    smiles = 'c1ccccc1'
    print(f'Output: {level_function(smiles)}')
