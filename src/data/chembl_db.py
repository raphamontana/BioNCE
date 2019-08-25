#   Script to gather data from Chembl Database using webresource_client
#   Resources and examples from the client usage can be found on https://github.com/chembl/chembl_webresource_client

from chembl_webresource_client.new_client import new_client
import pandas as pd

chembl_id = "CHEMBL3837"
std_type = "IC50"

class Chembl:
    """
    Target is Cathepsin-L (CHEMBL3837) and Standard Type is IC50
    Stores data into panda dataframe
    """
    @classmethod
    def wrc(cls):
        """
        wrc = web_resource_client
        Gather data from chembl_database using the wrc
        Returns a dataframe containing molecule chembl_id and its canonical smiles
        """
        chemblwrc = new_client.activity\
            .filter(target_chembl_id=chembl_id)\
            .filter(standard_type=std_type)\
            .only(['molecule_chembl_id','canonical_smiles'])

        return pd.DataFrame(chemblwrc)


#a = Chembl.wrc()
#print(a.shape)