# Standardize SMILES
import numpy as np
from chembl_structure_pipeline import standardizer
from rdkit import Chem

def standardize_smiles(smiles):
    try:
        # Convert to RDKit mol
        m = Chem.MolFromSmiles(smiles)
        # Neutralize and separate from counterion
        m = standardizer.get_parent_mol(m, neutralize=True, check_exclusion=True, verbose=False)[0]
        # Standardize representation
        m = standardizer.standardize_mol(m, check_exclusion=True)
        # Return RDKit canonical SMILES
        # remove stereochemistry information to avoid some activity cliffs
        return Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(m, isomericSmiles=False)))
    except:
        return np.nan
