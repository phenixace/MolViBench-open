from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol):
    """
    生成分子的 3D 构象。
    """
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        # 添加氢原子
        mol_h = Chem.AddHs(mol_obj)
        # 生成 3D 构象
        result = AllChem.EmbedMolecule(mol_h, AllChem.ETKDG())
        if result == -1:
            return None
        # 获取 3D 坐标
        conf = mol_h.GetConformer()
        coords = []
        for i in range(mol_h.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            coords.append((mol_h.GetAtomWithIdx(i).GetSymbol(), pos.x, pos.y, pos.z))
        return coords
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CCO"  # 乙醇
    result = level_function(smiles)
    if result:
        for atom_symbol, x, y, z in result:
            print(f"{atom_symbol}: ({x:.4f}, {y:.4f}, {z:.4f})")
