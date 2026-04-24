from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol1, mol2):
    """模拟芳香化合物的 Friedel-Crafts 酰基化。"""
    try:
        reaction_smarts = '[c:1][H].[C:2](=[O:3])[Cl]>>[c:1][C:2](=[O:3])'
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        aromatic = Chem.MolFromSmiles(mol1)
        acyl_halide = Chem.MolFromSmiles(mol2)
        if aromatic is None or acyl_halide is None:
            return None
        products = rxn.RunReactants((aromatic, acyl_halide))
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
    acyl_halide = "CC(=O)Cl"
    print(f"Friedel-Crafts 酰基化产物: {level_function(aromatic, acyl_halide)}")
