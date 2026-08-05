from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


PHARMACOPHORE_PATTERNS = {
    "Hydrogen-bond donor": "[N,O,S;!H0]",
    "Hydrogen-bond acceptor": "[N,O,S;!+]",
    "Negative charge center": "[$([O-,S-,P-]),$([OH]-[C,S,P]=O)]",
    "Positive charge center": "[+1,$([N;H2,H3])]",
    "Five-membered aromatic ring": "a1aaaa1",
    "Six-membered aromatic ring": "a1aaaaa1",
    "Hydrophobic center": "[C;!$(C-[!#6])]",
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


if __name__ == '__main__':
    smiles = 'CC(=O)Oc1ccccc1C(=O)O'
    result = level_function(smiles)
    if result:
        for p in result:
            print(f"Output: {p['type']}{p['count']}")
