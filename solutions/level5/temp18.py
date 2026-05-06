from rdkit import Chem

TOXIC_SMARTS = [
]

def level_function(mols):
    try:
        patterns = []
        for smarts in TOXIC_SMARTS:
            pat = Chem.MolFromSmarts(smarts)
            if pat:
                patterns.append(pat)

        results = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            is_toxic = False
            for pat in patterns:
                if mol.HasSubstructMatch(pat):
                    is_toxic = True
                    break
            if not is_toxic:
                results.append(Chem.MolToSmiles(mol))
        return results
    except Exception as e:
        print(e)
        return None
