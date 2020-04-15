from chembl_webresource_client.utils import utils
from chembl_webresource_client.new_client import new_client
from rdkit import DataStructs
from rdkit.Chem import PandasTools, RDKFingerprint
import pandas as pd
import numpy as np
from rdkit.Chem.PandasTools import ChangeMoleculeRendering
from rdkit.Chem.Fingerprints import FingerprintMols


class ChemblData:
    """"docstring """

    def __init__( self, name=None, smiles=None, chemblid=None, pchembl=None, activity=None, fps=None, sim=None ):
        self.name = name
        self.smiles = smiles
        self.chemblid = chemblid
        self.pchembl = pchembl
        self.activity = activity
        self.fps = fps
        self.sim = sim

    def get_data( self ):
        ids2df = pd.DataFrame()
        ids = ['CHEMBL3837', 'CHEMBL4072', 'CHEMBL3563', 'CHEMBL2954', 'CHEMBL268']
        for i in range( 0, len( ids ) ):
            activities = new_client.activity.filter( target_chembl_id__in=ids[i], standard_relation='=',
                                                    standard_type__iregex='(IC50|Ki)', pchembl_value__isnull=False )
            ids2df = ids2df.append( pd.DataFrame( activities.only( ['molecule_chembl_id', 'canonical_smiles',
                                                                 'standard_type', 'pchembl_value'] ) ) )
            print( '{}'.format( i ) )

        # some filtering that could only be applied after the request
        ids2df = ids2df.drop( columns=['value', 'type'] )
        ids2df['pchembl_value'] = ids2df['pchembl_value'].astype( float )
        ids2df = ids2df[ids2df.pchembl_value >= 4]
        ids2df['activity_class'] = np.where( ids2df['pchembl_value'] < 6, 'A',
                                            np.where( ids2df['pchembl_value'] < 8, 'B',
                                                     np.where( ids2df['pchembl_value'] <= 10, 'C', np.nan ) ) )
        ChangeMoleculeRendering( renderer='String' )
        ids2df = ids2df.dropna( subset=['canonical_smiles'] )
        #add fps column

        return ids2df

        #

a = ChemblData()
data = a.get_data()
print( data.shape )

