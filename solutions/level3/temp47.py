from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

TARGET_FEATURES = {
    "Kinase inhibitor": ["c1ccncc1", "c1ccc2[nH]ccc2c1", "[NH]C(=O)"],
    "GPCR ligand": ["c1ccc2c(c1)CCNC2", "[NH]c1ccccc1", "c1ccncc1"],
    "Protease inhibitor": ["C(=O)N", "[OH]C(=O)", "NC(=O)C"],
    "Nuclear receptor ligand": ["C1CCC2C(C1)CCC1C2CCC2(C)C1CCC2O", "c1ccc(O)cc1"],
    "Ion channel modulator": ["c1ccccc1N", "[NH2]CCCC", "C(=O)N"],
}

def level_function(mol):
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
