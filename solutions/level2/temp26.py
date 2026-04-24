from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol):
    """计算构象 RMSD。"""
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=2, randomSeed=42)
        if len(cids) < 2:
            return None
        for cid in cids:
            AllChem.MMFFOptimizeMolecule(mol, confId=cid)
        rmsd = AllChem.GetConformerRMS(mol, cids[0], cids[1])
        return rmsd
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccccc1CCO"
    result = level_function(smiles)
    print(f"两个构象的 RMSD: {result}")
