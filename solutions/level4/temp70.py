from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, FilterCatalog



BRENK_SMARTS = {
    "aldehyde": "[CH1](=O)", "epoxide": "C1OC1",
    "peroxide": "OO", "azide": "N=[N+]=[N-]",
    "disulfide": "SS", "nitro": "[N+](=O)[O-]",
}


def level_function(smiles_list):

    try:
        categories = {"safe": [], "risky": [], "reject": []}

        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue

            canonical = Chem.MolToSmiles(mol)
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
            num_atoms = mol.GetNumAtoms()
            mr = Crippen.MolMR(mol)


            veber = rot_bonds <= 10 and tpsa <= 140


            ghose = (160 <= mw <= 480 and -0.4 <= logp <= 5.6 and
                     40 <= num_atoms <= 480 and 20 <= mr <= 130)


            brenk_pass = True
            for name, smarts in BRENK_SMARTS.items():
                pattern = Chem.MolFromSmarts(smarts)
                if pattern and mol.HasSubstructMatch(pattern):
                    brenk_pass = False
                    break

            passes = [veber, ghose, brenk_pass]
            num_pass = sum(passes)

            if num_pass == 3:
                categories["safe"].append(canonical)
            elif num_pass == 0:
                categories["reject"].append(canonical)
            else:
                categories["risky"].append(canonical)

        return {
            "safe_count": len(categories["safe"]),
            "risky_count": len(categories["risky"]),
            "reject_count": len(categories["reject"]),
            "safe": categories["safe"],
            "risky": categories["risky"],
            "reject": categories["reject"]
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    mols = ['c1ccc(NC(=O)c2ccccc2)cc1', 'CCO', 'c1ccc([N+](=O)[O-])cc1', 'CC(C)Cc1ccc(C(C)C(=O)O)cc1', 'CCCCCCCCCCCCCCCCCCCC']
    result = level_function(mols)
    if result:
        print(f"Output: {result['safe_count']}{result['risky_count']}{result['reject_count']}")
