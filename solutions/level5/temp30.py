from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors


def level_function(mol):
    """给定一个候选分子，探索不同的芳环取代方式并比较 LogP。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol_obj)
        if num_aromatic_rings == 0:
            return None

        orig_logp = Descriptors.MolLogP(mol_obj)
        orig_smi = Chem.MolToSmiles(mol_obj)

        substituents = [
            ('[cH:1]>>[c:1]C', '甲基'),
            ('[cH:1]>>[c:1]O', '羟基'),
            ('[cH:1]>>[c:1]N', '氨基'),
            ('[cH:1]>>[c:1]F', '氟'),
            ('[cH:1]>>[c:1]Cl', '氯'),
            ('[cH:1]>>[c:1]OC', '甲氧基'),
            ('[cH:1]>>[c:1]C(F)(F)F', '三氟甲基'),
            ('[cH:1]>>[c:1][N+](=O)[O-]', '硝基'),
        ]

        results = []
        for smarts, name in substituents:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((mol_obj,))
            if products:
                product = products[0][0]
                try:
                    Chem.SanitizeMol(product)
                    smi = Chem.MolToSmiles(product)
                    logp = Descriptors.MolLogP(product)
                    results.append({
                        'substituent': name,
                        'smiles': smi,
                        'logp': round(logp, 2),
                        'delta_logp': round(logp - orig_logp, 2)
                    })
                except Exception:
                    continue

        results.sort(key=lambda x: x['logp'])
        return {"original_logp": round(orig_logp, 2), "derivatives": results}
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccccc1"
    result = level_function(smiles)
    if result:
        print(f"原始 LogP: {result['original_logp']}")
        for d in result['derivatives']:
            print(f"  {d['substituent']}: LogP={d['logp']} (Δ={d['delta_logp']})")
