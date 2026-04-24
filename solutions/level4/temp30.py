from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):
    """给定分子 → 判断是否含腈基 → 若有 → 水解为羧酸 → 计算 LogP。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含腈基
        pattern = Chem.MolFromSmarts('[C:1]#N')
        has_nitrile = mol_obj.HasSubstructMatch(pattern)

        if not has_nitrile:
            return None

        # Step 2: 水解腈基为羧酸
        rxn = AllChem.ReactionFromSmarts('[C:1]#N>>[C:1](=O)O')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        # Step 3: 计算 LogP
        logp = rdMolDescriptors.CalcCrippenDescriptors(product)[0]

        return {
            "has_nitrile": has_nitrile,
            "product": product_smiles,
            "logp": round(logp, 4)
        }
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CC#N"
    print(f"result: {level_function(smiles)}")
