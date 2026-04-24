from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol):
    """模拟醇的脱水生成烯烃。"""
    try:
        reaction_smarts = '[C:1]([OH])[C:2]([H])>>[C:1]=[C:2]'
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        molecule = Chem.MolFromSmiles(mol)
        if molecule is None:
            return None
        products = rxn.RunReactants((molecule,))
        if not products:
            return None
        result_smiles = []
        for product_set in products:
            for product in product_set:
                Chem.SanitizeMol(product)
                smi = Chem.MolToSmiles(product)
                if smi not in result_smiles:
                    result_smiles.append(smi)
        return result_smiles
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "CCO"
    print(f"醇脱水生成烯烃产物: {level_function(smiles)}")
