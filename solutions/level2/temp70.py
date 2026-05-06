from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.ML.Cluster import Butina

def level_function(smiles_list, distance_threshold=0.4):
    try:
        mols = []
        fps = []
        valid_smiles = []
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                mols.append(mol)
                fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048))
                valid_smiles.append(Chem.MolToSmiles(mol))

        if len(fps) < 2:
            return None

        n = len(fps)
        dists = []
        for i in range(1, n):
            for j in range(i):
                dist = 1 - DataStructs.TanimotoSimilarity(fps[i], fps[j])
                dists.append(dist)

        clusters = Butina.ClusterData(dists, n, distance_threshold, isDistData=True)

        mol_cluster = {}
        for cluster_id, cluster in enumerate(clusters):
            for mol_idx in cluster:
                mol_cluster[mol_idx] = cluster_id

        result = []
        for i, smi in enumerate(valid_smiles):
            result.append({
                "smiles": smi,
                "cluster_id": mol_cluster.get(i, -1)
            })

        return result
    except Exception as e:
        print(e)
        return None
