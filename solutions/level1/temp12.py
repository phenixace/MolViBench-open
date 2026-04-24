from rdkit import Chem

def level_function(mol):
    """
    输出分子的 canonical SMILES。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CCO"
    print(f"Canonical SMILES: {level_function(smiles)}")