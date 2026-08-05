from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold


def level_function(mol1, mol2):

    try:
        m1 = Chem.MolFromSmiles(mol1)
        m2 = Chem.MolFromSmiles(mol2)
        if m1 is None or m2 is None:
            return None


        scaf1 = MurckoScaffold.GetScaffoldForMol(m1)
        scaf2 = MurckoScaffold.GetScaffoldForMol(m2)
        scaf1_smi = Chem.MolToSmiles(scaf1)
        scaf2_smi = Chem.MolToSmiles(scaf2)



        hybrid1_mols = AllChem.ReplaceSubstructs(m2, scaf2, scaf1)

        hybrid2_mols = AllChem.ReplaceSubstructs(m1, scaf1, scaf2)

        hybrid1_smi = None
        hybrid2_smi = None

        if hybrid1_mols:
            try:
                Chem.SanitizeMol(hybrid1_mols[0])
                hybrid1_smi = Chem.MolToSmiles(hybrid1_mols[0])
            except Exception:
                pass

        if hybrid2_mols:
            try:
                Chem.SanitizeMol(hybrid2_mols[0])
                hybrid2_smi = Chem.MolToSmiles(hybrid2_mols[0])
            except Exception:
                pass


        fp1 = AllChem.GetMorganFingerprintAsBitVect(m1, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(m2, 2, nBits=2048)

        result = {
            "scaffold1": scaf1_smi,
            "scaffold2": scaf2_smi,
        }

        if hybrid1_smi:
            h1_mol = Chem.MolFromSmiles(hybrid1_smi)
            if h1_mol:
                h1_fp = AllChem.GetMorganFingerprintAsBitVect(h1_mol, 2, nBits=2048)
                result["hybrid1"] = {
                    "smiles": hybrid1_smi,
                    "sim_to_mol1": round(DataStructs.TanimotoSimilarity(fp1, h1_fp), 4),
                    "sim_to_mol2": round(DataStructs.TanimotoSimilarity(fp2, h1_fp), 4)
                }

        if hybrid2_smi:
            h2_mol = Chem.MolFromSmiles(hybrid2_smi)
            if h2_mol:
                h2_fp = AllChem.GetMorganFingerprintAsBitVect(h2_mol, 2, nBits=2048)
                result["hybrid2"] = {
                    "smiles": hybrid2_smi,
                    "sim_to_mol1": round(DataStructs.TanimotoSimilarity(fp1, h2_fp), 4),
                    "sim_to_mol2": round(DataStructs.TanimotoSimilarity(fp2, h2_fp), 4)
                }

        return result
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smi1 = 'c1ccc(NC(=O)C)cc1'
    smi2 = 'c1ccnc(O)c1'
    result = level_function(smi1, smi2)
    if result:
        print(f"Output: {result['scaffold1']}")
        print(f"Output: {result['scaffold2']}")
