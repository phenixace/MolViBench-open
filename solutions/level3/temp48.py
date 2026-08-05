from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


PBT_TOXIC_SMARTS = [
    "[$(a[N+](=O)[O-])]",
    "c1cc([Cl])c([Cl])c([Cl])c1",
    "[Sn]",
    "[Pb]",
    "[Hg]",
    "[Cd]",
]


def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        issues = []

        logp = Descriptors.MolLogP(mol_obj)
        if logp > 4.5:
            issues.append(f"LogP={round(logp, 2)} > 4.5, may be bioaccumulative")

        mw = Descriptors.MolWt(mol_obj)

        for smarts in PBT_TOXIC_SMARTS:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol_obj.HasSubstructMatch(pattern):
                issues.append(f"Contains a PBT-related toxic substructure: {smarts}")

        halogen_count = sum(1 for atom in mol_obj.GetAtoms()
                           if atom.GetAtomicNum() in [9, 17, 35, 53])
        if halogen_count >= 3:
            issues.append(f"Contains {halogen_count} halogen atoms and may be persistent")

        return {
            "compliant": len(issues) == 0,
            "issues": issues,
            "logp": round(logp, 2),
            "mw": round(mw, 2)
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1cc(Cl)c(Cl)c(Cl)c1'
    print(f'Output: {level_function(smiles)}')
