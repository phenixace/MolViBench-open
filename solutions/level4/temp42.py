from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


def level_function(mol):
    """给定分子 → 判断是否含炔烃 → 若有 → 卤化 → 计算分子式。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含炔烃
        pattern = Chem.MolFromSmarts('[C:1]#[C:2]')
        has_alkyne = mol_obj.HasSubstructMatch(pattern)

        if not has_alkyne:
            return None

        # Step 2: 卤化 (炔烃加溴)
        rxn = AllChem.ReactionFromSmarts('[C:1]#[C:2]>>[C:1](Br)=[C:2]Br')
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        # Step 3: 计算分子式
        formula = rdMolDescriptors.CalcMolFormula(product)

        return {
            "has_alkyne": has_alkyne,
            "product": product_smiles,
            "formula": formula
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "CC#C"
    print(f"result: {level_function(smiles)}")
