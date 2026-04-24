from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol):
    """
    给定苯，生成所有二取代甲基衍生物。
    """
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        # 先生成一取代甲基衍生物
        rxn = AllChem.ReactionFromSmarts('[cH:1]>>[c:1]C')
        first_products = rxn.RunReactants((mol_obj,))
        first_unique = set()
        first_mols = []
        for product_tuple in first_products:
            for product in product_tuple:
                Chem.SanitizeMol(product)
                smi = Chem.MolToSmiles(product)
                if smi not in first_unique:
                    first_unique.add(smi)
                    first_mols.append(product)
        # 再对每个一取代产物进行第二次取代
        unique_smiles = set()
        for m in first_mols:
            second_products = rxn.RunReactants((m,))
            for product_tuple in second_products:
                for product in product_tuple:
                    Chem.SanitizeMol(product)
                    smi = Chem.MolToSmiles(product)
                    unique_smiles.add(smi)
        return sorted(list(unique_smiles))
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "c1ccccc1"  # 苯
    result = level_function(smiles)
    print(f"二取代甲基衍生物: {result}")
