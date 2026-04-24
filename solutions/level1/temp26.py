from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

def level_function(mol):
    """
    判断分子是否符合 Lipinski Rule of Five。
    """
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
         # Lipinski Rule of Five 条件
         # 分子量 ≤ 500
         # LogP ≤ 5
         # 氢键供体数 ≤ 5
         # 氢键受体数 ≤ 10
        return (mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10)
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "CC[C@H](F)C(=O)O"
    print(f"是否符合 Lipinski Rule of Five: {level_function(smiles)}")