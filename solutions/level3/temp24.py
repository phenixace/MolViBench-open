from rdkit import Chem
from collections import Counter


def level_function(reactants, products):
    """判断反应是否平衡。"""
    try:
        def count_atoms(smiles_list):
            total = Counter()
            for smi in smiles_list:
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    return None
                mol = Chem.AddHs(mol)
                for atom in mol.GetAtoms():
                    total[atom.GetSymbol()] += 1
            return total

        reactant_atoms = count_atoms(reactants)
        product_atoms = count_atoms(products)
        if reactant_atoms is None or product_atoms is None:
            return None
        return reactant_atoms == product_atoms
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    reactants = ["CC(=O)O", "CCO"]
    products = ["CC(=O)OCC", "O"]
    print(f"反应是否平衡: {level_function(reactants, products)}")
