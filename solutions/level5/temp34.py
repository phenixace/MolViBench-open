from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def level_function(mols):
    try:
        mol_data = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            logp = Descriptors.MolLogP(mol)
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            mol_data.append({
                'smiles': Chem.MolToSmiles(mol),
                'logp': round(logp, 2),
                'tpsa': round(tpsa, 2),
            })

        pareto_front = []
        for i, d in enumerate(mol_data):
            dominated = False
            for j, d2 in enumerate(mol_data):
                if i == j:
                    continue
                if (d2['logp_dist'] <= d['logp_dist'] and d2['tpsa'] <= d['tpsa'] and
                        (d2['logp_dist'] < d['logp_dist'] or d2['tpsa'] < d['tpsa'])):
                    dominated = True
                    break
            if not dominated and d['tpsa'] < 120:
                pareto_front.append(d)

        pareto_front.sort(key=lambda x: x['logp_dist'])
        return pareto_front
    except Exception as e:
        print(e)
        return None
