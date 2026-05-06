from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

PHARMACOPHORE_PATTERNS = {
    "H-bond donor": "[#7H,#8H,#16H]",
    "H-bond acceptor": "[#7,#8,#16;!H0;v2,v3,v4,v5]",
    "Positive charge center": "[+,NH3+,NH2+,NH+]",
    "Negative charge center": "[-,C(=O)[O-],S(=O)(=O)[O-]]",
    "Aromatic ring": "a1aaaaa1",
    "Hydrophobic center": "[CH2,CH3,c]",
    "Halogen": "[F,Cl,Br,I]",
}

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        pharmacophores = []
        for name, smarts in PHARMACOPHORE_PATTERNS.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                continue
            matches = mol_obj.GetSubstructMatches(pattern)
            if matches:
                pharmacophores.append({
                    "type": name,
                    "count": len(matches),
                    "atom_indices": [list(m) for m in matches]
                })

        return pharmacophores if pharmacophores else None
    except Exception as e:
        print(e)
        return None
