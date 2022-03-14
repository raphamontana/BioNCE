# PAINS filter - RDKit
# https://projects.volkamerlab.org/teachopencadd/talktorials/T003_compound_unwanted_substructures.html
# https://www.rdkit.org/docs/source/rdkit.Chem.rdfiltercatalog.html

import pandas as pd
from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem import Draw
# from rdkit.Chem.Draw import IPythonConsole

# Read list with SMARTS for PAINS
pains_list = pd.read_csv("../pains/wehi_pains.csv", header=None)
pains_list.columns = ["SMARTS", "ID"]

# initialize filter
params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
catalog = FilterCatalog(params)


def pass_pains(smiles):
    # Returns True if mol has no PAINS substructures
    mol = Chem.MolFromSmiles(smiles)
    entry = catalog.GetFirstMatch(mol)  # Get first matching PAINS
    # entry.GetDescription().capitalize()
    if entry:
        return False
    return True


def pains_highlight(smiles):
    # Returns image of molecule with PAINS highlighted
    # Note: maybe some SMARTS strings will throw errors: https://github.com/rdkit/rdkit/issues/3912
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if not pass_pains(smiles):  # verifies if mol has PAINS
        match_list = []
        # Draw mol grid highlighting PAINS
        for smarts in pains_list["SMARTS"]:
            match = mol.GetSubstructMatch(Chem.MolFromSmarts(smarts))
            if match:
                match_list.append(match)
        ms = [mol]*len(match_list)
        img = Draw.MolsToGridImage(ms,
                                   molsPerRow=2,
                                   subImgSize=(400, 300),
                                   highlightAtomLists=match_list)
    return img
