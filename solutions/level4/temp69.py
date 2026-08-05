from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import sys
import os


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None


        rxn_templates = [
            '[cH:1]>>[c:1]C',
            '[cH:1]>>[c:1]O',
            '[cH:1]>>[c:1]N',
            '[cH:1]>>[c:1]F',
            '[cH:1]>>[c:1]Cl',
            '[cH:1]>>[c:1]OC',
            '[cH:1]>>[c:1]C(=O)O',
        ]

        derivatives = set()
        for rxn_sma in rxn_templates:
            rxn = AllChem.ReactionFromSmarts(rxn_sma)
            products = rxn.RunReactants((mol_obj,))
            for prod_set in products:
                for prod in prod_set:
                    try:
                        Chem.SanitizeMol(prod)
                        derivatives.add(Chem.MolToSmiles(prod))
                    except Exception:
                        pass


        lipinski_pass = []
        for smi in derivatives:
            m = Chem.MolFromSmiles(smi)
            if m is None:
                continue
            mw = Descriptors.MolWt(m)
            logp = Descriptors.MolLogP(m)
            hbd = Descriptors.NumHDonors(m)
            hba = Descriptors.NumHAcceptors(m)
            if mw < 500 and logp < 5 and hbd <= 5 and hba <= 10:
                lipinski_pass.append((smi, m))

        if not lipinski_pass:
            return {"derivatives_total": len(derivatives), "lipinski_pass": 0, "top3": []}


        try:
            from rdkit.Chem import RDConfig
            sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
            import sascorer
            sa_func = sascorer.calculateScore
        except Exception:

            sa_func = lambda m: Descriptors.BertzCT(m) / 100.0

        scored = []
        for smi, m in lipinski_pass:
            sa = sa_func(m)
            scored.append({"smiles": smi, "SA_Score": round(sa, 4)})

        scored.sort(key=lambda x: x["SA_Score"])

        return {
            "derivatives_total": len(derivatives),
            "lipinski_pass": len(lipinski_pass),
            "top3": scored[:3]
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccccc1'
    result = level_function(smiles)
    if result:
        print(f"Output: {result['derivatives_total']}{result['lipinski_pass']}")
        for r in result['top3']:
            print(f"Output: {r['smiles']}{r['SA_Score']}")
