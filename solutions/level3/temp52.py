from rdkit import Chem

def level_function(reaction_smiles):
    try:
        parts = reaction_smiles.split(">>")
        if len(parts) != 2:
            return None

        reactant_smis = parts[0].split(".")
        product_smis = parts[1].split(".")

        reactants = [Chem.MolFromSmiles(s) for s in reactant_smis if Chem.MolFromSmiles(s)]
        products = [Chem.MolFromSmiles(s) for s in product_smis if Chem.MolFromSmiles(s)]

        if not reactants or not products:
            return None

        r_heavy = sum(m.GetNumHeavyAtoms() for m in reactants)
        p_heavy = sum(m.GetNumHeavyAtoms() for m in products)

        r_bonds = sum(m.GetNumBonds() for m in reactants)
        p_bonds = sum(m.GetNumBonds() for m in products)

        r_rings = sum(Chem.rdMolDescriptors.CalcNumRings(m) for m in reactants)
        p_rings = sum(Chem.rdMolDescriptors.CalcNumRings(m) for m in products)

        halide_pattern = Chem.MolFromSmarts("[F,Cl,Br,I]")
        r_has_halide = any(m.HasSubstructMatch(halide_pattern) for m in reactants)
        p_has_halide = any(m.HasSubstructMatch(halide_pattern) for m in products)

        db_pattern = Chem.MolFromSmarts("[#6]=[#6]")
        r_has_db = any(m.HasSubstructMatch(db_pattern) for m in reactants)
        p_has_db = any(m.HasSubstructMatch(db_pattern) for m in products)

        if r_heavy == p_heavy and len(reactants) == 1 and len(products) == 1:
            return "Rearrangement"
        elif r_has_halide and not p_has_halide and len(products) >= len(reactants):
            return "Nucleophilic Substitution"
        elif not r_has_db and p_has_db and p_bonds < r_bonds:
            return "Elimination"
        elif r_has_db and not p_has_db:
            return "Electrophilic Addition"
        elif r_has_halide and not p_has_halide:
            return "Nucleophilic Substitution"
        elif len(reactants) > len(products):
            return "Addition"
        elif len(reactants) < len(products):
            return "Elimination"
        else:
            return "Substitution"
    except Exception as e:
        print(e)
        return None
