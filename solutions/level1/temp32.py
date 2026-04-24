from rdkit import Chem

def level_function(mol):
    """
    将分子转为原子特征矩阵。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        atom_features = []
        for atom in mol.GetAtoms():
            features = [
                atom.GetAtomicNum(),  # 原子序数
                atom.GetDegree(),     # 原子度
                atom.GetFormalCharge(), # 形式电荷
                atom.GetHybridization().real, # 杂化类型
                int(atom.GetIsAromatic()) # 是否芳香族
            ]
            atom_features.append(features)
        return atom_features
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CC[C@H](F)C(=O)O"
    print(f"原子特征矩阵: {level_function(smiles)}")