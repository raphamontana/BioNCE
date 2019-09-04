# Import statements
import os
import pandas as pd

# rdkit:
# from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import RDKFingerprint
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw, PandasTools
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols


class PrepareData(object):
    """

    """

    def read_from_csv(self):
        """
        Read .csv on /data/raw/ folder
        This returns a dataframe for every .csv contained in the raw folder
        """
        # set path to read from .csv downloaded files
        path_to_data = os.path.normpath(os.getcwd() + os.sep + 'raw/')
        csv_targets = os.listdir(path_to_data)
        csv_targets = filter(lambda x: x.endswith(".csv"), csv_targets)
        df_returned = []
        for element in csv_targets:
            df = pd.read_csv('' + path_to_data + '/{}'.format(element), sep=',', encoding='utf-8', index_col='id')
            df.name = '{}'.format(str(element))
            df_returned.append(df)
        return df_returned

    def get_fps(self, df):
        """" This will only work if there is no missing values on smilesCol.
        Add a column named ROMol (by default) with the molecule object containing its fingerprint"""
        PandasTools.AddMoleculeColumnToFrame(df, smilesCol='canonical_smiles', includeFingerprints=True)

    def update_csv(self, df):
        path_to_download = os.path.normpath(os.getcwd() + os.sep + '/raw/')
        df.to_csv(r'' + path_to_download + '/{}.csv'.format(str(df.name)),
                  sep=',', encoding='utf-8', index=True, index_label='id')

    def drop_nan_smiles(self, df):
        df = df.dropna(subset=['canonical_smiles'])
        return df


def main():
    data = PrepareData()
    df1, df2 = data.read_from_csv()
    df1 = data.drop_nan_smiles(df1)
    df2 = data.drop_nan_smiles(df2)
    data.get_fps(df1)
    data.get_fps(df2)
    print(df1.shape, df2.shape)
    df1
    df2

if __name__ == '__main__':
    main()
