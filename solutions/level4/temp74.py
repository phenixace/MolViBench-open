from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, FilterCatalog

BRENK_SMARTS = {
    "aldehyde": "[CH1](=O)", "epoxide": "C1OC1", "peroxide": "OO",
    "azide": "N=[N+]=[N-]", "disulfide": "SS", "nitro": "[N+](=O)[O-]",
}

def level_function(smiles_list):
    try:
        mols = []
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                mols.append((Chem.MolToSmiles(mol), mol))

        stages = [{"stage": "input", "count": len(mols)}]

        lipinski_pass = []
        for smi, mol in mols:
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            hbd = rdMolDescriptors.CalcNumHBD(mol)
            hba = rdMolDescriptors.CalcNumHBA(mol)
            if mw < 500 and logp < 5 and hbd <= 5 and hba <= 10:
                lipinski_pass.append((smi, mol))
        stages.append({"stage": "Lipinski", "count": len(lipinski_pass)})

        veber_pass = []
        for smi, mol in lipinski_pass:
            rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
            tpsa = Descriptors.TPSA(mol)
            if rot <= 10 and tpsa <= 140:
                veber_pass.append((smi, mol))
        stages.append({"stage": "Veber", "count": len(veber_pass)})

        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog.FilterCatalog(params)

        pains_pass = []
        for smi, mol in veber_pass:
            if catalog.GetFirstMatch(mol) is None:
                pains_pass.append((smi, mol))
        stages.append({"stage": "PAINS", "count": len(pains_pass)})

        brenk_pass = []
        for smi, mol in pains_pass:
            clean = True
            for name, smarts in BRENK_SMARTS.items():
                pat = Chem.MolFromSmarts(smarts)
                if pat and mol.HasSubstructMatch(pat):
                    clean = False
                    break
            if clean:
                brenk_pass.append((smi, mol))
        stages.append({"stage": "Brenk", "count": len(brenk_pass)})

        return {
            "cascade_results": stages,
            "surviving_smiles": [smi for smi, _ in brenk_pass]
        }
    except Exception as e:
        print(e)
        return None
