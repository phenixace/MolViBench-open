from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def level_function(mols):
    """给定一组分子，筛选旋转键数 < 10 的分子。"""
    try:
        results = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
            if rot < 10:
                results.append({
                    "smiles": Chem.MolToSmiles(mol),
                    "rotatable_bonds": rot
                })
        return results
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = ["CCO", "c1ccccc1", "CCCCCCCCCCCCCCCCCCC"]
    print(f"旋转键 < 10: {level_function(smiles_list)}")
