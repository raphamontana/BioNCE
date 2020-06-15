from chembl_webresource_client.new_client import new_client
import pandas as pd
import numpy as np
import os
import sys

# Set path to download .csv files to "raw" project folder
path_to_download = os.path.normpath( os.getcwd() + os.sep + os.pardir + os.sep + os.pardir + '/data/raw/' )


class ChemblDataExtraction( object ):
    """CDE stands for Chembl Data Extraction.
    The purpose of this class is to use the chembls official webclient to
    gather and filter data based on certain targets and standard types
    """

    def acquire( self, chembl_id ):
        """ chembl_id for primary targets: Cruzipain('CHEMBL3563') and Cathepsin L ('CHEMBL3837')
        standard types expected: 'IC50', 'Ki'
        The dataframe will filter all molecules that has <4 pchembl value and also classify its intervals based on:
            Columns: 'pchembl_interval' A ={4:5,9}; B={6:7,9}; C={8:10}
        """

        chemblwrc = new_client.activity \
            .filter( target_chembl_id=chembl_id ) \
            .filter( standard_relation='=' ) \
            .filter( standard_type__iregex='(IC50|Ki)' ) \
            .only( ['molecule_chembl_id', 'canonical_smiles', 'standard_type', 'pchembl_value'] )

        # Create dataframe from the web query above
        chemblwrc = pd.DataFrame( chemblwrc )

        chemblwrc = chemblwrc.drop( columns=['value', 'type'] )

        # Convert p_chembl value type to float so that it can be filtered by a threshold
        chemblwrc['pchembl_value'] = chemblwrc['pchembl_value'].astype( float )
        chemblwrc = chemblwrc[chemblwrc.pchembl_value >= 4]

        # create new column and assign a label based on p_chembl interval
        chemblwrc['interval'] = np.where( chemblwrc['pchembl_value'] < 6, 'A',
                                         np.where( chemblwrc['pchembl_value'] < 8, 'B',
                                                  np.where( chemblwrc['pchembl_value'] <= 10, 'C', np.nan ) ) )
        chemblwrc.name = chembl_id
        return( chemblwrc )

    def write( self, dataframe ):
        dataframe.to_csv( r'' + path_to_download + '/{}.csv'.format( str( dataframe.name ) ), sep=',', encoding='utf-8',
                         index=True, index_label='id' )

