from rdkit import Chem

def level_function(mol1, mol2):
    """
    判断两个 SMILES 是否等价。
    """
    try:
        mol1 = Chem.MolFromSmiles(mol1)
        if mol1 is None:
            return None
        mol2 = Chem.MolFromSmiles(mol2)
        if mol2 is None:
            return None
        return Chem.MolToSmiles(mol1) == Chem.MolToSmiles(mol2)
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles1 = "CC[C@H](F)C(=O)O"
    smiles2 = "CC[C@@H](F)C(=O)O"
    print(f"是否等价: {level_function(smiles1, smiles2)}")