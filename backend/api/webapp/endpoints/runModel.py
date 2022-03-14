"""
BioNCE API Script.

Implement a Flask API to calculate molecules descriptors with RDKIT and run a
Machine Learning model.
"""

###############################################################################
# Import libs.
###############################################################################

import joblib
import pandas as pd
from flask import Flask, jsonify, request
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


###############################################################################
# Load ML models.
###############################################################################

models = dict()
models["catb_model-v22.2"] = joblib.load("../models/catb_model-v22.2.pkl")
models["catl_model-v22.2"] = joblib.load("../models/catl_model-v22.2.pkl")
models["cats_model-v22.2"] = joblib.load("../models/cats_model-v22.2.pkl")
models["mpro_model-v22.2"] = joblib.load("../models/mpro_model-v22.2.pkl")


###############################################################################
# Calculate descriptors function.
###############################################################################

def _calculateDescriptors(mol):
    df = pd.DataFrame(index=[0])
    df["SlogP"] = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    df["SMR"] = rdMolDescriptors.CalcCrippenDescriptors(mol)[1]
    df["LabuteASA"] = rdMolDescriptors.CalcLabuteASA(mol)
    df["TPSA"] = Descriptors.TPSA(mol)
    df["AMW"] = Descriptors.MolWt(mol)
    df["ExactMW"] = rdMolDescriptors.CalcExactMolWt(mol)
    df["NumLipinskiHBA"] = rdMolDescriptors.CalcNumLipinskiHBA(mol)
    df["NumLipinskiHBD"] = rdMolDescriptors.CalcNumLipinskiHBD(mol)
    df["NumRotatableBonds"] = rdMolDescriptors.CalcNumRotatableBonds(mol)
    df["NumHBD"] = rdMolDescriptors.CalcNumHBD(mol)
    df["NumHBA"] = rdMolDescriptors.CalcNumHBA(mol)
    df["NumAmideBonds"] = rdMolDescriptors.CalcNumAmideBonds(mol)
    df["NumHeteroAtoms"] = rdMolDescriptors.CalcNumHeteroatoms(mol)
    df["NumHeavyAtoms"] = Chem.rdchem.Mol.GetNumHeavyAtoms(mol)
    df["NumAtoms"] = Chem.rdchem.Mol.GetNumAtoms(mol)
    df["NumRings"] = rdMolDescriptors.CalcNumRings(mol)
    df["NumAromaticRings"] = rdMolDescriptors.CalcNumAromaticRings(mol)
    df["NumSaturatedRings"] = rdMolDescriptors.CalcNumSaturatedRings(mol)
    df["NumAliphaticRings"] = rdMolDescriptors.CalcNumAliphaticRings(mol)
    df["NumAromaticHeterocycles"] = \
        rdMolDescriptors.CalcNumAromaticHeterocycles(mol)
    df["NumSaturatedHeterocycles"] = \
        rdMolDescriptors.CalcNumSaturatedHeterocycles(mol)
    df["NumAliphaticHeterocycles"] = \
        rdMolDescriptors.CalcNumAliphaticHeterocycles(mol)
    df["NumAromaticCarbocycles"] = \
        rdMolDescriptors.CalcNumAromaticCarbocycles(mol)
    df["NumSaturatedCarbocycles"] = \
        rdMolDescriptors.CalcNumSaturatedCarbocycles(mol)
    df["NumAliphaticCarbocycles"] = \
        rdMolDescriptors.CalcNumAliphaticCarbocycles(mol)
    df["FractionCSP3"] = rdMolDescriptors.CalcFractionCSP3(mol)
    df["Chi0v"] = rdMolDescriptors.CalcChi0v(mol)
    df["Chi1v"] = rdMolDescriptors.CalcChi1v(mol)
    df["Chi2v"] = rdMolDescriptors.CalcChi2v(mol)
    df["Chi3v"] = rdMolDescriptors.CalcChi3v(mol)
    df["Chi4v"] = rdMolDescriptors.CalcChi4v(mol)
    df["Chi1n"] = rdMolDescriptors.CalcChi1n(mol)
    df["Chi2n"] = rdMolDescriptors.CalcChi2n(mol)
    df["Chi3n"] = rdMolDescriptors.CalcChi3n(mol)
    df["Chi4n"] = rdMolDescriptors.CalcChi4n(mol)
    df["HallKierAlpha"] = rdMolDescriptors.CalcHallKierAlpha(mol)
    df["kappa1"] = rdMolDescriptors.CalcKappa1(mol)
    df["kappa2"] = rdMolDescriptors.CalcKappa2(mol)
    df["kappa3"] = rdMolDescriptors.CalcKappa3(mol)
    slogp_VSA = list(map(lambda i: "slogp_VSA"+str(i), list(range(1, 13))))
    df = df.assign(**dict(zip(slogp_VSA, rdMolDescriptors.SlogP_VSA_(mol))))
    smr_VSA = list(map(lambda i: "smr_VSA"+str(i), list(range(1, 11))))
    df = df.assign(**dict(zip(smr_VSA, rdMolDescriptors.SMR_VSA_(mol))))
    peoe_VSA = list(map(lambda i: "peoe_VSA"+str(i), list(range(1, 15))))
    df = df.assign(**dict(zip(peoe_VSA, rdMolDescriptors.PEOE_VSA_(mol))))
    MQNs = list(map(lambda i: "MQN"+str(i), list(range(1, 43))))
    df = df.assign(**dict(zip(MQNs, rdMolDescriptors.MQNs_(mol))))
    return df


###############################################################################
# Run prediction function.
###############################################################################

def prediction(model, smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return("SMILES Parse Error")
    df = _calculateDescriptors(mol)
    if model in models:
        res = models[model].predict(df)
        return res[0]
    else:
        return("Model not found.")
