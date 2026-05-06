from rdkit import Chem

def level_function(mol):
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        smiles_set = set()
        for _ in range(100):
            smi = Chem.MolToRandomSmilesVect(mol, 1)[0]
            smiles_set.add(smi)
        smiles_set.add(Chem.MolToSmiles(mol))
        return list(smiles_set)
    except Exception as e:
        print(e)
        return None
