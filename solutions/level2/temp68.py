from rdkit import Chem

AA_SMILES = {
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
