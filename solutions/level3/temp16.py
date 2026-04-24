from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol1, mol2):
    """模拟芳香化合物的 Friedel-Crafts 烷基化。"""
    try:
        reaction_smarts = '[c:1][H].[C:2][Cl]>>[c:1][C:2]'
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        aromatic = Chem.MolFromSmiles(mol1)
        alkyl_halide = Chem.MolFromSmiles(mol2)
        if aromatic is None or alkyl_halide is None:
            return None
        products = rxn.RunReactants((aromatic, alkyl_halide))
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
    aromatic = "c1ccccc1"
    alkyl_halide = "CCl"
    print(f"Friedel-Crafts 烷基化产物: {level_function(aromatic, alkyl_halide)}")
