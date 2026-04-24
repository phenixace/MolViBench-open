from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(reactants, reaction_smarts):
    """判断反应是否可能发生。"""
    try:
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        mols = []
        for smi in reactants:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                return False
            mols.append(mol)
        products = rxn.RunReactants(tuple(mols))
        return len(products) > 0
    except Exception as e:
        print(e)
        return False


if __name__ == "__main__":
    reactants = ["CC(=O)O", "CCN"]
    reaction_smarts = "[C:1](=O)[OH].[N:2]>>[C:1](=O)[N:2]"
    print(f"反应是否可能发生: {level_function(reactants, reaction_smarts)}")
