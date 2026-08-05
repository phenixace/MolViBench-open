from rdkit import Chem



AA_SMILES = {
    'A': 'C',
    'R': 'CCCCNC(=N)N',
    'N': 'CC(=O)N',
    'D': 'CC(=O)O',
    'C': 'CS',
    'E': 'CCC(=O)O',
    'Q': 'CCC(=O)N',
    'G': '[H]',
    'H': 'Cc1cnc[nH]1',
    'I': 'C(CC)C',
    'L': 'CC(C)C',
    'K': 'CCCCN',
    'M': 'CCSC',
    'F': 'Cc1ccccc1',
    'P': '',
    'S': 'CO',
    'T': 'C(O)C',
    'W': 'Cc1c[nH]c2ccccc12',
    'Y': 'Cc1ccc(O)cc1',
    'V': 'C(C)C',
}


def level_function(sequence):

    try:
        if not sequence or not isinstance(sequence, str):
            return None

        sequence = sequence.upper().strip()







        residues = []
        for aa in sequence:
            if aa not in AA_SMILES:
                return None
            residues.append(aa)

        if not residues:
            return None


        parts = []
        for i, aa in enumerate(residues):
            side_chain = AA_SMILES[aa]
            if aa == 'G':
                parts.append("NCC(=O)")
            elif aa == 'P':

                parts.append("N1CCCC1C(=O)")
            else:
                parts.append(f"N[C@@H]({side_chain})C(=O)")


        peptide_smiles = "".join(parts) + "O"


        mol = Chem.MolFromSmiles(peptide_smiles)
        if mol is None:

            peptide_smiles_no_stereo = peptide_smiles.replace("[C@@H]", "C")
            mol = Chem.MolFromSmiles(peptide_smiles_no_stereo)
            if mol is None:
                return None
            return Chem.MolToSmiles(mol)

        return Chem.MolToSmiles(mol)
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    seq = 'AGL'
    print(f'Output: {seq}{level_function(seq)}')
    seq2 = 'ACDEF'
    print(f'Output: {seq2}{level_function(seq2)}')
