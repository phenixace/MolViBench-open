from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions


def level_function(mol):
    """给定分子，生成所有手性中心构型。"""
    try:
        molecule = Chem.MolFromSmiles(mol)
        if molecule is None:
            return None

        opts = StereoEnumerationOptions(tryEmbedding=True, unique=True, onlyUnassigned=False)
        isomers = list(EnumerateStereoisomers(molecule, options=opts))

        result_smiles = []
        for isomer in isomers:
            smi = Chem.MolToSmiles(isomer, isomericSmiles=True)
            if smi not in result_smiles:
                result_smiles.append(smi)

        return result_smiles if result_smiles else [Chem.MolToSmiles(molecule, isomericSmiles=True)]
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    mol = "CC(O)C(F)Cl"
    print(f"所有手性中心构型: {level_function(mol)}")
