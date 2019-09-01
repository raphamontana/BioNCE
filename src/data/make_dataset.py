from load_chembl_data import ChemblDataExtraction
import pandas as pd

def main():
    """ Write .csv file for each target containing both IC50 and Ki types
    """
    data = ChemblDataExtraction()
    cathepsin = data.acquire('CHEMBL3837')
    cruzipain = data.acquire('CHEMBL3563')
    data.write(cathepsin)
    data.write(cruzipain)
if __name__ == '__main__':
    main()
