"""
Module to verify if tested compound is within the applicability domain of the model

Modules: numpy, pandas, pickle, rdkit
Files: dataframe with limits for the molecular descriptors, dataframes of training sets (fingerprints)

descriptor_ad: returns True if at least 5 descriptors are within the set limits

fingerprint_ad: returns True if the fingerprint for the tested compound has a Tanimoto score >= 0.5 with at least 5% of the dataset
"""

import numpy as np
import pandas as pd
import pickle
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

# Import dataframe with descriptor limits
descs_df = pickle.load(open("../ad/descs_ad.pkl", "rb"))

# Import fingerprint lists
fps_catb = pickle.load(open("../ad/fps_catb.pkl", "rb"))
fps_catl = pickle.load(open("../ad/fps_catl.pkl", "rb"))
fps_cats = pickle.load(open("../ad/fps_cats.pkl", "rb"))
fps_mpro = pickle.load(open("../ad/fps_mpro.pkl", "rb"))


def descriptor_ad(smiles, enzyme):
    # input: SMILES of compound to be tested and any of these strings: "catb", "catl", "cats", "mpro"
    mol = Chem.MolFromSmiles(smiles)
    descs = ["MolLogP", "MolWt", "NumHAcceptors",
             "NumHDonors", "NumRotatableBonds", "TPSA"]
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(descs)
    mol_descs = calc.CalcDescriptors(mol)
    mol_descs_dict = dict(zip(descs, mol_descs))

    # Compare to dataset
    fail = 0
    for desc in descs:
        if (mol_descs_dict[desc] < descs_df.loc[enzyme, desc][0]) or (mol_descs_dict[desc] > descs_df.loc[enzyme, desc][1]):
            fail += 1

    return fail < 1


def fingerprint_ad(smiles, enzyme_fp_list, score=0.17):
    # input: SMILES of compound to be tested and list of fingerprints as RDKit ExplicitBitVect objects
    mol = Chem.MolFromSmiles(smiles)
    fp_mol = Chem.AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    sim = []

    for fp in enzyme_fp_list:
        sim.append(DataStructs.FingerprintSimilarity(
            fp_mol, fp, metric=DataStructs.TanimotoSimilarity))

    sim_pass = np.array(sim) >= score
    if (sim_pass.sum() > len(sim)/20):
        return True
    return False
