from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem, Draw


def level_function(mol1, mol2):
    """基于 MCS 对两个分子进行 2D 对齐并可视化。"""
    try:
        m1 = Chem.MolFromSmiles(mol1)
        m2 = Chem.MolFromSmiles(mol2)
        if m1 is None or m2 is None:
            return None

        # Find MCS
        mcs_result = rdFMCS.FindMCS([m1, m2], timeout=10)
        mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
        if mcs_mol is None:
            return None

        # Generate 2D coords
        AllChem.Compute2DCoords(m1)
        AllChem.Compute2DCoords(m2)

        # Align m2 to m1 based on MCS
        match1 = m1.GetSubstructMatch(mcs_mol)
        match2 = m2.GetSubstructMatch(mcs_mol)

        if not match1 or not match2:
            return None

        AllChem.GenerateDepictionMatching2DStructure(m2, m1,
            atomMap=list(zip(match2, match1)))

        # Generate SVG
        svg = Draw.MolsToGridImage([m1, m2], molsPerRow=2,
                                    subImgSize=(400, 300),
                                    legends=["Mol1", "Mol2"],
                                    useSVG=True)
        return {
            "mcs_smarts": mcs_result.smartsString,
            "mcs_numAtoms": mcs_result.numAtoms,
            "svg": svg
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smi1 = "c1ccc(CCO)cc1"
    smi2 = "c1ccc(CCN)cc1"
    result = level_function(smi1, smi2)
    if result:
        print(f"MCS SMARTS: {result['mcs_smarts']}, Atoms: {result['mcs_numAtoms']}")
