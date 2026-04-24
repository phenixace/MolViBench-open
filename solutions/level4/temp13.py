from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs

def level_function(mol):
    """给定分子 → 判断是否含手性中心 → 若有 → 生成对映体 → 计算相似度。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含手性中心
        Chem.AssignStereochemistry(mol_obj, cleanIt=True, force=True)
        chiral_centers = Chem.FindMolChiralCenters(mol_obj)
        has_chiral = len(chiral_centers) > 0

        if not has_chiral:
            return None

        # Step 2: 生成对映体（翻转所有手性中心）
        enantiomer = Chem.RWMol(mol_obj)
        for atom in enantiomer.GetAtoms():
            chiral_tag = atom.GetChiralTag()
            if chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
            elif chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)

        enantiomer_mol = enantiomer.GetMol()
        enantiomer_smiles = Chem.MolToSmiles(enantiomer_mol)

        # Step 3: 计算 Tanimoto 相似度
        fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2)
        fp1 = fpgen.GetFingerprint(mol_obj)
        fp2 = fpgen.GetFingerprint(enantiomer_mol)
        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)

        return {
            "has_chiral": has_chiral,
            "product": enantiomer_smiles,
            "similarity": round(similarity, 4)
        }
    except Exception as e:
        print(e)
        return None

if __name__ == "__main__":
    smiles = "C[C@H](O)CC"
    print(f"result: {level_function(smiles)}")
