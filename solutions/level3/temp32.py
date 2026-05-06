from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

def level_function(mol):
    try:
        molecule = Chem.MolFromSmiles(mol)
        if molecule is None:
            return None

        opts = StereoEnumerationOptions(tryEmbedding=True, unique=True)
        isomers = list(EnumerateStereoisomers(molecule, options=opts))

        result_smiles = []
        for isomer in isomers:
            smi = Chem.MolToSmiles(isomer, isomericSmiles=True)
            if smi not in result_smiles:
                result_smiles.append(smi)

        return result_smiles if result_smiles else [Chem.MolToSmiles(molecule)]
    except Exception as e:
        print(e)
        return None
