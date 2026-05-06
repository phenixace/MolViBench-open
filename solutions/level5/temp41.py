from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem.Scaffolds import MurckoScaffold

def level_function(mols):
    try:
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
