from rdkit import Chem

def level_function(mol):
    """
    判断分子是否是手性分子。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        return len(chiral_centers) > 0
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CC[C@H](F)C(=O)O"
    print(f"是否是手性分子: {level_function(smiles)}")