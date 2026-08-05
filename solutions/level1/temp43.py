from rdkit import Chem

FUNCTIONAL_GROUPS = {
    "(-OH)": "[OX2H]",
    "(-NH2)": "[NX3;H2]",
    "(-COOH)": "C(=O)[OX2H1]",
    "(-CHO)": "[CX3H1](=O)[#6]",
    "(C=O)": "[CX3](=O)[#6]",
    "(-COOR)": "C(=O)O[#6]",
    "(F/Cl/Br/I)": "[F,Cl,Br,I]",
    "(Ar)": "a",
    "(P)": "[PX4]",
    "(-SO2-)": "S(=O)(=O)[#6]"
}

def level_function(mol):
    try:
        mol = Chem.MolFromSmiles(mol) if isinstance(mol, str) else mol
        if mol is None:
            return None

        found_groups = []
        for name, smarts in FUNCTIONAL_GROUPS.items():
            patt = Chem.MolFromSmarts(smarts)
            if mol.HasSubstructMatch(patt):
                found_groups.append(name)

        return found_groups if found_groups else None

    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "PC(N)C[C@H](F)C(=O)O"
    print(f"Output: {level_function(smiles)}")
