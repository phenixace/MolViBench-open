from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors, Lipinski

def level_function(mols):
    try:
        mol_data = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Lipinski.NumHDonors(mol)
            hba = Lipinski.NumHAcceptors(mol)
            passes_lipinski = mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10
            mol_data.append({
                'smiles': Chem.MolToSmiles(mol),
                'fp': fp,
                'passes_lipinski': passes_lipinski
            })

        if len(mol_data) < 2:
            return None

        results = []
        for i, d in enumerate(mol_data):
            if not d['passes_lipinski']:
                continue
            is_diverse = True
            for j, d2 in enumerate(mol_data):
                if i == j:
                    continue
                sim = DataStructs.TanimotoSimilarity(d['fp'], d2['fp'])
                if sim >= 0.4:
                    is_diverse = False
                    break

            if is_diverse:
                results.append(d['smiles'])

        if not results:
            for d in mol_data:
                if d['passes_lipinski']:
                    avg_sim = sum(
                        DataStructs.TanimotoSimilarity(d['fp'], d2['fp'])
                        for d2 in mol_data if d2['smiles'] != d['smiles']
                    ) / (len(mol_data) - 1)
                    if avg_sim < 0.4:
                        results.append(d['smiles'])

        return results
    except Exception as e:
        print(e)
        return None
