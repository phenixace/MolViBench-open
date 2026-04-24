from rdkit import Chem


# Standard amino acid SMILES (L-amino acids, neutral form)
AA_SMILES = {
    'A': 'C',            # Alanine
    'R': 'CCCCNC(=N)N',  # Arginine
    'N': 'CC(=O)N',      # Asparagine
    'D': 'CC(=O)O',      # Aspartic acid
    'C': 'CS',           # Cysteine
    'E': 'CCC(=O)O',     # Glutamic acid
    'Q': 'CCC(=O)N',     # Glutamine
    'G': '[H]',          # Glycine (no side chain)
    'H': 'Cc1cnc[nH]1',  # Histidine
    'I': 'C(CC)C',       # Isoleucine
    'L': 'CC(C)C',       # Leucine
    'K': 'CCCCN',        # Lysine
    'M': 'CCSC',         # Methionine
    'F': 'Cc1ccccc1',    # Phenylalanine
    'P': '',             # Proline (special - cyclic)
    'S': 'CO',           # Serine
    'T': 'C(O)C',        # Threonine
    'W': 'Cc1c[nH]c2ccccc12',  # Tryptophan
    'Y': 'Cc1ccc(O)cc1',      # Tyrosine
    'V': 'C(C)C',        # Valine
}


def level_function(sequence):
    """将氨基酸单字母序列转换为对应的线性肽 SMILES（L-氨基酸，中性形式）。"""
    try:
        if not sequence or not isinstance(sequence, str):
            return None

        sequence = sequence.upper().strip()

        # Build peptide using RDKit
        # Construct the linear peptide SMILES manually
        # Each amino acid: N[C@@H](side_chain)C(=O)
        # Peptide bond: -NH-CH(R)-CO-NH-CH(R)-CO-...
        # N-terminus: free NH2, C-terminus: free COOH

        residues = []
        for aa in sequence:
            if aa not in AA_SMILES:
                return None
            residues.append(aa)

        if not residues:
            return None

        # Build peptide SMILES
        parts = []
        for i, aa in enumerate(residues):
            side_chain = AA_SMILES[aa]
            if aa == 'G':
                parts.append("NCC(=O)")
            elif aa == 'P':
                # Proline is special (cyclic side chain)
                parts.append("N1CCCC1C(=O)")
            else:
                parts.append(f"N[C@@H]({side_chain})C(=O)")

        # Join with peptide bonds and add termini
        peptide_smiles = "".join(parts) + "O"

        # Validate
        mol = Chem.MolFromSmiles(peptide_smiles)
        if mol is None:
            # Try without stereochemistry
            peptide_smiles_no_stereo = peptide_smiles.replace("[C@@H]", "C")
            mol = Chem.MolFromSmiles(peptide_smiles_no_stereo)
            if mol is None:
                return None
            return Chem.MolToSmiles(mol)

        return Chem.MolToSmiles(mol)
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    seq = "AGL"
    print(f"序列 {seq} -> SMILES: {level_function(seq)}")
    seq2 = "ACDEF"
    print(f"序列 {seq2} -> SMILES: {level_function(seq2)}")
