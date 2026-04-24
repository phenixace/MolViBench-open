from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol):
    """生成多个构象并筛选 RMSD 差异 >0.5 的。"""
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=50, params=AllChem.ETKDG())
        if len(conf_ids) == 0:
            return None
        selected = [conf_ids[0]]
        for i in range(1, len(conf_ids)):
            dominated = False
            for j in selected:
                rmsd = AllChem.GetConformerRMS(mol, j, conf_ids[i])
                if rmsd <= 0.5:
                    dominated = True
                    break
            if not dominated:
                selected.append(conf_ids[i])
        return selected
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccccc1CCO"
    result = level_function(smiles)
    print(f"RMSD > 0.5 的构象 ID: {result}")
