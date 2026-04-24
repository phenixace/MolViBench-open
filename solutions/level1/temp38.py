from rdkit import Chem

def level_function(mol, substructure):
    """
    给定 SMARTS 模式，判断是否匹配。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        substructure = Chem.MolFromSmarts(substructure)
        
        if substructure is None:
            return None
        
        return mol.HasSubstructMatch(substructure)
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "C(N)C[C@H](F)C(=O)O"
    print(f"给定 SMARTS 模式，判断是否匹配。: {level_function(smiles, "[NX3;H2]")}")