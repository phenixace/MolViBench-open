from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


# 常见靶点的关键子结构特征 SMARTS
TARGET_FEATURES = {
    "激酶抑制剂": ["c1ccncc1", "c1ccc2[nH]ccc2c1", "[NH]C(=O)"],
    "GPCR 配体": ["c1ccc2c(c1)CCNC2", "[NH]c1ccccc1", "c1ccncc1"],
    "蛋白酶抑制剂": ["C(=O)N", "[OH]C(=O)", "NC(=O)C"],
    "核受体配体": ["C1CCC2C(C1)CCC1C2CCC2(C)C1CCC2O", "c1ccc(O)cc1"],
    "离子通道调节剂": ["c1ccccc1N", "[NH2]CCCC", "C(=O)N"],
}


def level_function(mol):
    """给定分子，预测可能的靶点活性。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        predictions = []
        for target, smarts_list in TARGET_FEATURES.items():
            match_count = 0
            for smarts in smarts_list:
                pattern = Chem.MolFromSmarts(smarts)
                if pattern and mol_obj.HasSubstructMatch(pattern):
                    match_count += 1
            if match_count > 0:
                confidence = round(match_count / len(smarts_list), 2)
                predictions.append({
                    "target": target,
                    "matched_features": match_count,
                    "total_features": len(smarts_list),
                    "confidence": confidence
                })

        predictions.sort(key=lambda x: x["confidence"], reverse=True)
        return predictions if predictions else None
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccc2c(c1)cc(CC(=O)O)[nH]2"
    result = level_function(smiles)
    if result:
        for p in result:
            print(f"  {p['target']}: 置信度 {p['confidence']}")
