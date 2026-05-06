from rdkit import Chem

FUNCTIONAL_GROUPS = {
    "Hydroxyl (-OH)": "[OX2H]", 
    "Amino (-NH2)": "[NX3;H2]", 
    "Carboxyl (-COOH)": "C(=O)[OX2H1]", 
    "Aldehyde (-CHO)": "[CX3H1](=O)[#6]", 
    "Ketone (C=O)": "[CX3](=O)[#6]", 
    "Ester (-COOR)": "C(=O)O[#6]", 
    "Halogen (F/Cl/Br/I)": "[F,Cl,Br,I]", 
    "Aromatic ring (Ar)": "a", 
    "Phosphine (P)": "[PX4]", 
    "Sulfonyl (-SO2-)": "S(=O)(=O)[#6]"
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
