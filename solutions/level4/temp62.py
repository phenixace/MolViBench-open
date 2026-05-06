from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdFMCS

def level_function(smiles_list):
    try:
        mols = []
        fps = []
        valid_smiles = []
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                mols.append(mol)
                fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048))
                valid_smiles.append(Chem.MolToSmiles(mol))

        if len(mols) < 2:
            return None

        n = len(mols)
        best_sim = -1
        best_pair = (0, 1)
        sim_matrix = []
        for i in range(n):
            row = []
            for j in range(n):
                if i == j:
                    row.append(1.0)
                elif j > i:
                    sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
                    row.append(round(sim, 4))
                    if sim > best_sim:
                        best_sim = sim
                        best_pair = (i, j)
                else:
                    row.append(sim_matrix[j][i])
            sim_matrix.append(row)

        i, j = best_pair
        mcs_result = rdFMCS.FindMCS([mols[i], mols[j]], timeout=10)

        return {
            "num_molecules": n,
            "most_similar_pair": {
                "mol1": valid_smiles[i],
                "mol2": valid_smiles[j],
                "similarity": round(best_sim, 4)
            },
            "MCS": {
                "smarts": mcs_result.smartsString,
                "numAtoms": mcs_result.numAtoms,
                "numBonds": mcs_result.numBonds
            },
            "similarity_matrix": sim_matrix
        }
    except Exception as e:
        print(e)
        return None
