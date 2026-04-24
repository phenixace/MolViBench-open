from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol):
    """给定分子 → 判断是否含醛基 → 若有 → 氧化为羧酸 → 计算 TPSA。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含醛基
        pattern = Chem.MolFromSmarts('[CX3H1](=O)')
        has_aldehyde = mol_obj.HasSubstructMatch(pattern)

        if not has_aldehyde:
            return None

        # Step 2: 氧化醛基为羧酸
        rxn = AllChem.ReactionFromSmarts('[CH:1]=O>>[C:1](=O)O')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        # Step 3: 计算 TPSA
        tpsa = rdMolDescriptors.CalcTPSA(product)

        return {
            "has_aldehyde": has_aldehyde,
            "product": product_smiles,
            "tpsa": round(tpsa, 4)
        }
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CC=O"
    print(f"result: {level_function(smiles)}")
