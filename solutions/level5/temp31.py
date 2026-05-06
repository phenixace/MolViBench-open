from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors

def level_function(mols):
    try:
        mol_data = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            qed = Descriptors.qed(mol)
            mol_data.append({
                'smiles': Chem.MolToSmiles(mol),
                'fp': fp,
                'qed': round(qed, 4)
            })

        if len(mol_data) < 2:
            return None

        results = []
        for i, d in enumerate(mol_data):
            sims = []
            for j, d2 in enumerate(mol_data):
                if i != j:
                    sim = DataStructs.TanimotoSimilarity(d['fp'], d2['fp'])
                    sims.append(sim)
            avg_sim = sum(sims) / len(sims) if sims else 0
            results.append({
                'smiles': d['smiles'],
                'avg_similarity': round(avg_sim, 4),
                'qed': d['qed']
            })

        return results
    except Exception as e:
        print(e)
        return None
