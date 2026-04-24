from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol):
    """模拟苯的硝化反应。"""
    try:
        reaction_smarts = '[c:1][H]>>[c:1][N+](=O)[O-]'
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
    smiles = "c1ccccc1"
    print(f"硝化反应产物: {level_function(smiles)}")
