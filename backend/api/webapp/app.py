"""
BioNCE API Script.

Implement a Flask API to calculate molecules descriptors with RDKIT and run a
Machine Learning model.
"""

###############################################################################
# Import libs.
###############################################################################

from functions.runModel import prediction
from functions.ad import fingerprint_ad, descriptor_ad
from functions.descriptors import descriptors
from functions.exmol_apply import get_similar_actives
from functions.filter_Ro5 import rule_of_five
from functions.pains import pass_pains, pains_highlight
from functions.standardize import standardize_smiles
from functions.fingerprint import getFingerprint

from flask import Flask, jsonify, request


###############################################################################
# Configure Flask app.
###############################################################################

app = Flask(__name__)


@app.route('/prediction')
def prediction_endpoint():
    """
    Prediction API.

    Returns
    -------
        TYPE: DESCRIPTION.

    """
    model = request.args.get('model')
    smiles = request.args.get('smiles')
    molPrediction = prediction(model, smiles)
    return jsonify({'prediction': repr(molPrediction)})


@app.route('/similar_actives')
def similar_actives_endpoint():
    base = request.args.get('base')
    model = request.args.get('model')
    similar_actives = get_similar_actives(base, model)
    return jsonify({'similar_actives': repr(similar_actives)})


@app.route('/standardize_smiles')
def standardize_smiles_endpoint():
    smiles = request.args.get('smiles')
    standardized = standardize_smiles(smiles)
    return jsonify({'standardize_smiles': repr(standardized)})


@app.route('/descriptor_ad')
def descriptor_ad_endpoint():
    smiles = request.args.get('smiles')
    enzyme = request.args.get('enzyme')
    molDescriptor_ad = descriptor_ad(smiles, enzyme)
    return jsonify({'descriptor_ad': repr(molDescriptor_ad)})


@app.route('/fingerprint_ad')
def fingerprint_ad_endpoint():
    smiles = request.args.get('smiles')
    enzyme_fp_list = request.args.get('enzyme_fp_list')
    fingerprints = fingerprint_ad(smiles, enzyme_fp_list)
    return jsonify({'fingerprint_ad': repr(fingerprints)})


@app.route('/descriptors')
def descriptors_endpoint():
    smiles = request.args.get('smiles')
    molDescriptors = descriptors(smiles)
    return jsonify({'descriptors': repr(molDescriptors)})


@app.route('/rule_of_five')
def rule_of_five_endpoint():
    values = request.args.get('values')
    ro5pass = rule_of_five(values)
    return jsonify({'rule_of_five': repr(ro5pass)})


@app.route('/pass_pains')
def pass_pains_endpoint():
    smiles = request.args.get('smiles')
    pains = pass_pains(smiles)
    return jsonify({'pass_pains': repr(pains)})


@app.route('/pains_highlight')
def pains_highlight_endpoint():
    smiles = request.args.get('smiles')
    pains = pains_highlight(smiles)
    return jsonify({'pains_highlight': repr(pains)})


@app.route('/fingerprint')
def fingerprint_endpoint():
    smiles = request.args.get('smiles')
    fingerprint = getFingerprint(smiles)
    return jsonify({'fingerprint': repr(fingerprint)})


if __name__ == '__main__':
    app.run(host='0.0.0.0')
