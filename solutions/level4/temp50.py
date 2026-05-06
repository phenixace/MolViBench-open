from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdMolDescriptors
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        pattern = Chem.MolFromSmarts('c1ccccc1')
        has_benzene = mol_obj.HasSubstructMatch(pattern)

        if not has_benzene:
            return None

        opts = StereoEnumerationOptions(unique=True)
        isomers = list(EnumerateStereoisomers(mol_obj, options=opts))

        if len(isomers) < 1:
            return None

        isomer_smiles = []
        fps = []
        for iso in isomers:
            smi = Chem.MolToSmiles(iso)
            if smi not in isomer_smiles:
                isomer_smiles.append(smi)
                fp = AllChem.GetMorganFingerprintAsBitVect(iso, 2, nBits=2048)
                fps.append(fp)

        n = len(fps)
        matrix = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i == j:
                    matrix[i][j] = 1.0
                elif j > i:
                    sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
                    matrix[i][j] = round(sim, 4)
                    matrix[j][i] = round(sim, 4)

        return {
            "has_benzene": has_benzene,
            "isomers": isomer_smiles,
            "similarity_matrix": matrix
        }
    except Exception as e:
        print(e)
        return None
