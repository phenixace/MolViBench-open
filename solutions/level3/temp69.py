from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors


def level_function(mol):
    """对分子同时应用多个过滤规则（Veber + Ghose + Egan）并输出每个规则的通过/未通过结果。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        mw = Descriptors.MolWt(mol_obj)
        logp = Crippen.MolLogP(mol_obj)
        tpsa = Descriptors.TPSA(mol_obj)
        rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol_obj)
        num_atoms = mol_obj.GetNumAtoms()
        mr = Crippen.MolMR(mol_obj)

        # Veber rules: RotBonds ≤ 10 AND TPSA ≤ 140
        veber = rot_bonds <= 10 and tpsa <= 140

        # Ghose rules: 160 ≤ MW ≤ 480, -0.4 ≤ LogP ≤ 5.6, 40 ≤ atoms ≤ 480, 20 ≤ MR ≤ 130
        ghose = (160 <= mw <= 480 and -0.4 <= logp <= 5.6 and
                 40 <= num_atoms <= 480 and 20 <= mr <= 130)

        # Egan rules: TPSA ≤ 131.6, LogP ≤ 5.88
        egan = tpsa <= 131.6 and logp <= 5.88

        return {
            "properties": {
                "MW": round(mw, 2),
                "LogP": round(logp, 2),
                "TPSA": round(tpsa, 2),
                "RotBonds": rot_bonds,
                "NumAtoms": num_atoms,
                "MR": round(mr, 2)
            },
            "Veber": veber,
            "Ghose": ghose,
            "Egan": egan,
            "passes_all": veber and ghose and egan
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccc(NC(=O)c2ccc(F)cc2)cc1"
    result = level_function(smiles)
    if result:
        print(f"Veber: {result['Veber']}, Ghose: {result['Ghose']}, Egan: {result['Egan']}")
        print(f"All pass: {result['passes_all']}")
