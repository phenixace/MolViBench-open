from rdkit import Chem

def level_function(mol):
    """
    输出分子中所有键类型。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return set([bond.GetBondType() for bond in mol.GetBonds()])
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CCOc1ccccc1"
    print(f"所有键类型: {level_function(smiles)}")