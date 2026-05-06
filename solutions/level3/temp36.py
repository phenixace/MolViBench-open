from rdkit import Chem
from rdkit.Chem import Descriptors

def level_function(mol1, mol2):
    try:
        m1 = Chem.MolFromSmiles(mol1)
        m2 = Chem.MolFromSmiles(mol2)
        if m1 is None or m2 is None:
            return False

        formula1 = Descriptors.MolecularFormula(m1)
        formula2 = Descriptors.MolecularFormula(m2)
        if formula1 != formula2:
            return False

        smi1_no_stereo = Chem.MolToSmiles(m1, isomericSmiles=False)
        smi2_no_stereo = Chem.MolToSmiles(m2, isomericSmiles=False)
        if smi1_no_stereo != smi2_no_stereo:
            return False

        smi1_stereo = Chem.MolToSmiles(m1, isomericSmiles=True)
        smi2_stereo = Chem.MolToSmiles(m2, isomericSmiles=True)
        if smi1_stereo == smi2_stereo:
            return False

        chiral1 = Chem.FindMolChiralCenters(m1, includeUnassigned=False)
        chiral2 = Chem.FindMolChiralCenters(m2, includeUnassigned=False)

        if len(chiral1) < 2 or len(chiral1) != len(chiral2):
            return False

        inverted_count = 0
        for (idx1, tag1), (idx2, tag2) in zip(sorted(chiral1), sorted(chiral2)):
            if idx1 != idx2:
                return False
            if tag1 != tag2:
                inverted_count += 1

        if inverted_count == len(chiral1):
            return False

        return inverted_count > 0
    except Exception as e:
        print(e)
        return False
