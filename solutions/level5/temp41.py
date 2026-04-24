from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem.Scaffolds import MurckoScaffold


def level_function(mols):
    """从一组已知药物中提取 scaffold → 生成取代衍生物 → 筛选 Lipinski 符合的分子。"""
    try:
        # Step 1: 提取 scaffold
        scaffolds = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            try:
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                scaffold_smi = Chem.MolToSmiles(scaffold)
                if scaffold_smi not in [s['smiles'] for s in scaffolds]:
                    scaffolds.append({'smiles': scaffold_smi, 'mol': scaffold})
            except Exception:
                continue

        if not scaffolds:
            return None

        # Step 2: 生成取代衍生物
        rxns = [
            '[cH:1]>>[c:1]C',
            '[cH:1]>>[c:1]O',
            '[cH:1]>>[c:1]F',
            '[cH:1]>>[c:1]Cl',
            '[cH:1]>>[c:1]N',
        ]

        derivatives = set()
        for scaffold in scaffolds:
            for rxn_smarts in rxns:
                rxn = AllChem.ReactionFromSmarts(rxn_smarts)
                products = rxn.RunReactants((scaffold['mol'],))
                for product_set in products:
                    for product in product_set:
                        try:
                            Chem.SanitizeMol(product)
                            smi = Chem.MolToSmiles(product)
                            derivatives.add(smi)
                        except Exception:
                            continue

        # Step 3: 筛选 Lipinski
        results = []
        for smi in derivatives:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Lipinski.NumHDonors(mol)
            hba = Lipinski.NumHAcceptors(mol)
            if mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10:
                results.append({
                    'smiles': smi,
                    'mw': round(mw, 2),
                    'logp': round(logp, 2)
                })

        return results[:20]
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = ["CC(=O)Oc1ccccc1C(=O)O", "c1ccc(O)cc1"]
    result = level_function(smiles_list)
    if result:
        for r in result[:5]:
            print(f"  {r['smiles']}: MW={r['mw']}, LogP={r['logp']}")
