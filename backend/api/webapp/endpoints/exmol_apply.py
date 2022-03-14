import exmol


def get_similar_actives(base, model):
    """
    base: a SMILES string
    model: a classification model that outputs 0 (inactive) / 1 (active)
    output: image from exmol
    """

    # Molecule of interest
    # base = 'Cc1onc(-c2ccccc2Cl)c1C(=O)NC1C(=O)N2C1SC(C)(C)C2C(=O)O'

    # Model - must output 0 for mol of interest
    def my_model(x):
        if x == base:
            return 0
        return model(x)

    sample_space = exmol.sample_space(base, lambda smi, sel: my_model(smi),
                                      batched=False, preset="medium")

    cfs5 = exmol.cf_explain(sample_space, nmols=5)

    return exmol.plot_cf(cfs5)
