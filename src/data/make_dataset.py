import load_chembl_data
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator


def main():
    """ Write .csv file for each target containing both IC50 and Ki types
    """
    data = load_chembl_data.ChemblDataExtraction()
    cathepsin = data.acquire('CHEMBL3837')
    cruzipain = data.acquire('CHEMBL3563')
    data.write(cathepsin)
    data.write(cruzipain)
    cruzipain
    cathepsin

main()
