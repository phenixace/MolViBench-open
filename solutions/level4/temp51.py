from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdMolDescriptors

def level_function(mol):
    try:
        from rdkit.Chem import Descriptors, BRICS, Crippen

        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        mw = Descriptors.MolWt(mol_obj)

        if mw > 300:
            frags = BRICS.BRICSDecompose(mol_obj)
            if not frags:
                return None
            best_frag = None
            best_size = 0
            for frag_smi in frags:
                frag_mol = Chem.MolFromSmiles(frag_smi)
                if frag_mol and frag_mol.GetNumHeavyAtoms() > best_size:
                    best_size = frag_mol.GetNumHeavyAtoms()
                    best_frag = frag_smi
            if best_frag is None:
                return None
            final_mol = Chem.MolFromSmiles(best_frag)
        else:
            rxn = AllChem.ReactionFromSmarts('[cH:1]>>[c:1]c1ccccc1')
            products = rxn.RunReactants((mol_obj,))
            if products:
                final_mol = products[0][0]
                Chem.SanitizeMol(final_mol)
            else:
                rxn2 = AllChem.ReactionFromSmarts('[CH:1]>>[C:1]c1ccccc1')
                products = rxn2.RunReactants((mol_obj,))
                if products:
                    final_mol = products[0][0]
                    Chem.SanitizeMol(final_mol)
                else:
                    final_mol = mol_obj

        final_smi = Chem.MolToSmiles(final_mol)
        logp = Crippen.MolLogP(final_mol)

        return {
            "original_MW": round(mw, 2),
            "action": "BRICS_decompose" if mw > 300 else "fragment_grow",
            "final_smiles": final_smi,
            "final_LogP": round(logp, 4)
        }
    except Exception as e:
        print(e)
        return None
