from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol):
    """模拟亲电取代反应。"""
    try:
        molecule = Chem.MolFromSmiles(mol)
        if molecule is None:
            return None
        reaction_smarts = '[cH:1]>>[c:1]Br'
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        products = rxn.RunReactants((molecule,))
        if not products:
            return None
        result_smiles = []
        for product_set in products:
            for product in product_set:
                try:
                    Chem.SanitizeMol(product)
                    smi = Chem.MolToSmiles(product)
                    if smi not in result_smiles:
                        result_smiles.append(smi)
                except Exception:
                    continue
        return result_smiles if result_smiles else None
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccccc1"
    print(f"亲电取代产物: {level_function(smiles)}")
