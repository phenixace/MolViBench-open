from rdkit import Chem
from rdkit.Chem import rdchem

def get_metal_atomic_numbers():
    """
    返回所有金属元素的原子序数集合。
    包括碱金属、碱土金属、过渡金属、镧系和锕系等。
    """
    # 手动定义金属范围（覆盖常见金属元素）
    metals = set()
    # 碱金属: Li, Na, K, Rb, Cs, Fr
    metals.update([3, 11, 19, 37, 55, 87])
    # 碱土金属: Be, Mg, Ca, Sr, Ba, Ra
    metals.update([4, 12, 20, 38, 56, 88])
    # 过渡金属: Sc (21) ~ Zn (30), Y (39) ~ Cd (48), Hf (72) ~ Hg (80)
    metals.update(range(21, 31))
    metals.update(range(39, 49))
    metals.update(range(72, 81))
    # 镧系: La (57) ~ Lu (71)
    metals.update(range(57, 72))
    # 锕系: Ac (89) ~ Lr (103)
    metals.update(range(89, 104))
    # 其他常见金属: Al (13), Ga (31), In (49), Tl (81), Sn (50), Pb (82), Bi (83), Po (84)
    metals.update([13, 31, 49, 50, 81, 82, 83, 84])
    return metals

def level_function(mol):
    """
    判断分子是否含有金属元素。
    """
    try:
        mol = Chem.MolFromSmiles(mol) if isinstance(mol, str) else mol
        if mol is None:
            return None

        metals = get_metal_atomic_numbers()
        metal_atoms = []

        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() in metals:
                metal_atoms.append(atom.GetSymbol())

        return True if metal_atoms else False

    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles1 = "PC(N)C[C@H](F)C(=O)O"   # 没有金属
    smiles2 = "[Na+].[O-]C(=O)C"       # 含Na金属
    smiles3 = "[Cu+2].[O-]C(=O)C"      # 含Cu金属
    print(f"{smiles1} 是否含有金属元素: {level_function(smiles1)}")
    print(f"{smiles2} 是否含有金属元素: {level_function(smiles2)}")
    print(f"{smiles3} 是否含有金属元素: {level_function(smiles3)}")
