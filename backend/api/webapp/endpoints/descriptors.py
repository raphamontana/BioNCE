# Calculate descriptors
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.ML.Descriptors import MoleculeDescriptors


def descriptors(smiles):
    """
    Calculate 8 most important descriptors
    input: SMILES (passed standardize_smiles())
    output: list with 8 descriptors calculated by RDKit, in order of desc_list
    """

    # SMILES must have passed standardize_smiles()
    mol = Chem.MolFromSmiles(smiles)

    # Descriptors list
    desc_list = ["MolLogP", "MolWt", "NumHAcceptors", "NumHDonors",
                 "NumRotatableBonds", "FractionCSP3", "TPSA", "NumAromaticRings"]
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(desc_list)
    values = calc.CalcDescriptors(mol)

    # Return list with values in the same order as declared in desc_list
    return values
