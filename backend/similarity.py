""" class definitions for preprocessing target data
"""
import os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import PandasTools, rdFingerprintGenerator
from rdkit.Chem.Fingerprints import FingerprintMols

def update_csv(df):
    path_to_download = os.path.normpath(os.getcwd() + os.sep + '/raw/')
    df.to_csv(r'' + path_to_download + '/{}'.format(str(df.name)),
              sep=',', encoding='utf-8', index=True, index_label='id')

class ChemblDataBatch(object):
    """  base class that reads the .csv files on /raw/ folder using os.paths
        Uses similiarity threshold = 1 for now
    """

    def __init__(self, data1=None, data2=None, threshold=1, **kwargs):
        self.sim_threshold = threshold
        path_to_data = os.path.normpath(os.getcwd() + os.sep + 'raw/')
        csv_targets = os.listdir(path_to_data)
        csv_targets = filter(lambda x: x.startswith("CHEMBL"), csv_targets)
        df_returned = []
        for element in csv_targets:
            df = pd.read_csv('' + path_to_data + '/{}'.format(element), sep=',', encoding='utf-8', index_col='id')
            # Droping null values from canonical_smiles
            df = df.dropna(subset=['canonical_smiles'])
            df.name = '{}'.format(str(element))
            df_returned.append(df)
            print('Read: {}'.format(element))

        self.number_of_elements = (len(df_returned))
        # TODO what if theres more than 2 here:
        self.data0 = df_returned[0]
        self.data1 = df_returned[1]
        self.sim_list = []


    def similarity_list(self, **kwargs):
        """ used to reset our internal state so that iteration
          starts again from the beginning
        """
        # Get dataframe names
        dataframes_read = []
        for x in range((int(self.number_of_elements))):
            dataframes_read.append('data{}'.format(x))
        # Iterate over multiple files
        data_to_compare = []
        c_smiles = []
        for a in dataframes_read:
            data_to_compare.append(getattr(self, a))
        for i, a in enumerate(data_to_compare):
            for ds in data_to_compare[i]['canonical_smiles']:
                cs = Chem.CanonSmiles(ds)
                c_smiles.append(cs)
            molecules = ([Chem.MolFromSmiles(x) for x in data_to_compare[i]['canonical_smiles']])
            fps = [FingerprintMols.FingerprintMol(x) for x in molecules]

            # the list for the dataframe
            qu, ta, sim = [], [], []

        # compare all fp pairwise without duplicates on the same csv
            for n in range(len(fps)):
                s = DataStructs.BulkTanimotoSimilarity(fps[n],
                                                       fps[n + 1:])  # +1 compare with the next to the last fp
                # collect the SMILES and values
                for m in range(len(s)):
                    if s[m] == self.sim_threshold:
                        qu.append(c_smiles[n])
                        ta.append(c_smiles[n + 1:][m])
                        sim.append(s[m])

            d = {'query_{}'.format(data_to_compare[i].name): qu,
                 'target_{}'.format(data_to_compare[i].name): ta,
                 'Similarity': sim}
            d = pd.DataFrame(data=d)
            d.name = 'sim_{}'.format(data_to_compare[i].name)
            self.sim_list.append(d)
            update_csv(d)

            # TODO make an abstract way of generating similarity for both files
            c_smiles1 = []
            c_smiles2 = []
            for ds in data_to_compare[0]['canonical_smiles']:
                cs = Chem.CanonSmiles(ds)
                c_smiles1.append(cs)
            for ds in data_to_compare[1]['canonical_smiles']:
                cs = Chem.CanonSmiles(ds)
                c_smiles2.append(cs)
            molecules1 = ([Chem.MolFromSmiles(x) for x in data_to_compare[0]['canonical_smiles']])
            molecules2 = ([Chem.MolFromSmiles(x) for x in data_to_compare[1]['canonical_smiles']])
            fps1 = [FingerprintMols.FingerprintMol(x) for x in molecules1]
            fps2 = [FingerprintMols.FingerprintMol(x) for x in molecules2]
            # the list for the dataframe
            qu, ta, sim = [], [], []

            # compare all fp pairwise without duplicates on the same csv
            for n in range(len(fps1)):
                s = DataStructs.BulkTanimotoSimilarity(fps1[n],
                                                       fps2[:])
                # collect the SMILES and values
                for m in range(len(s)):
                    if s[m] == self.sim_threshold:
                        qu.append(c_smiles1[n])
                        ta.append(c_smiles1[n])
                        sim.append(s[m])
            name_file_both_sim = 'both_files.csv'
            d = {'query_{}'.format(name_file_both_sim): qu,
                 'target_{}'.format(name_file_both_sim): ta,
                 'Similarity': sim}
            d = pd.DataFrame(data=d)
            d.name = 'sim_{}'.format(name_file_both_sim)
            self.sim_list.append(d)
            update_csv(d)


obj = ChemblDataBatch()
ChemblDataBatch.similarity_list(obj)
print((obj.sim_list[0].name))