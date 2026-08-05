from rdkit import Chem



TOXIC_SMARTS = [
    "[$(a[N+](=O)[O-]),$(a[N](=O)=O)]",
    "[NH2]a",
    "C1OC1",
    "[CX3](=[OX1])[F,Cl,Br,I]",
    "[N]=[N]",
    "OO",
    "[CX4]([F,Cl,Br,I])([F,Cl,Br,I])",
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


if __name__ == '__main__':
    smiles_list = ['CCO', 'c1ccc([N+](=O)[O-])cc1', 'c1ccccc1']
    print(f'Output: {level_function(smiles_list)}')
