from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def level_function(mols):
    try:
        scaffolds = set()
        valid_count = 0

        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            valid_count += 1
            try:
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                scaffold_smi = Chem.MolToSmiles(scaffold)
                scaffolds.add(scaffold_smi)
            except Exception:
                continue

        if valid_count == 0:
            return None

        diversity_index = len(scaffolds) / valid_count

        return {
            'num_molecules': valid_count,
            'num_unique_scaffolds': len(scaffolds),
            'scaffold_diversity_index': round(diversity_index, 4),
            'scaffolds': list(scaffolds)
        }
    except Exception as e:
        print(e)
        return None
