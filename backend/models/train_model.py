"""
IDS model training.

Montanari & Montanari
24/01/2020
"""

# Importing Modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle

from datetime import datetime

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

from sklearn import metrics
from sklearn.metrics import roc_curve
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import OneHotEncoder

"""
TODO Verificar necessidade destes imports
from rdkit import rdBase
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import PandasTools, RDKFingerprint
from rdkit.Chem.PandasTools import ChangeMoleculeRendering
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
help('modules')
"""


###############################################################################
# Data
processed = "../../data/processed/"
models = "../../data/models/"
file_in = processed + "meusdados.csv"
file_out = processed + "resultados.rdkit.csv"
meusdados = pd.read_csv(file_in)


###############################################################################
# Select all RDKit descriptors
calc = MoleculeDescriptors.MolecularDescriptorCalculator(
        [x[0] for x in Descriptors._descList]
    )
descriptorNames = calc.GetDescriptorNames()


###############################################################################
# Read smiles file and calculate the features for each molecule.
database = dict()
for a_mol in meusdados.values:
    chembl_id = a_mol[2]
    smiles = a_mol[1]
    m = Chem.MolFromSmiles(smiles)
    ds = calc.CalcDescriptors(m)
    database[chembl_id] = [smiles] + list(ds)
df = pd.DataFrame(database).T

###############################################################################
# Rename header and fill label column for the classifier.
col_names = ["SMILES"]
for name in descriptorNames:
    col_names.append(name)
df.columns = col_names

activity = []
for smiles in df.SMILES:
    filtro = meusdados[meusdados.SMILES == smiles]
    activity.append(filtro.LABEL.values[0])
df["LABEL"] = activity


###############################################################################
# Save training file.
df.to_csv(file_out, index=None)


###############################################################################
# Start training.
###############################################################################
df = pd.read_csv(file_out)                                   # Loading dataset.

""" TODO crossvalidate """

###############################################################################
# Preprocessing
# Prevalidate data.
df = df.replace(0, np.nan)                        # Remove columns with zeroes.
df = df.dropna(how='all', axis=0)
df = df.replace(np.nan, 0)
training_input = df.loc[:, descriptorNames]

###############################################################################
# Setting expected output.
training_output = []
for label in df.LABEL:
    if label == "Active":
        training_output.append(1)
    else:
        training_output.append(0)

###############################################################################
# Training step (10 fold validation).
print(datetime.now())
X_train, X_test, y_train, y_test = train_test_split(
    training_input, training_output, test_size=0.5, random_state=0
    )
clf = RandomForestClassifier(n_estimators=5000)              # Declaring Model.
gbClf = GradientBoostingClassifier()                         # Fitting Model.
scores = cross_val_score(clf, training_input, training_output, cv=10)

fpr, tpr, thresholds = metrics.roc_curve(training_output, scores, pos_label=2)

print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
print(datetime.now())

###############################################################################
# Save the model to disk.
pickle.dump(clf, open(models + "RandomForestClassifier.sav", "wb"))


###############################################################################
# Load the model from disk.
model = pickle.load(open(models + "RandomForestClassifier.sav", "rb"))
model.predict(training_input)

###############################################################################

""" TODO Incluir classificadores RF, SVM, MLP. """


# Supervised transformation based on gradient boosted trees
grd = GradientBoostingClassifier(n_estimators=10)
grd_enc = OneHotEncoder()
grd_lm = LogisticRegression(max_iter=1000)

grd.fit(X_train, y_train)
grd_enc.fit(grd.apply(X_train)[:, :, 0])
grd_lm.fit(grd_enc.transform(grd.apply(X_test)[:, :, 0]), y_test)

y_pred_grd_lm = grd_lm.predict_proba(
    grd_enc.transform(grd.apply(X_test)[:, :, 0]))[:, 1]
fpr_grd_lm, tpr_grd_lm, _ = roc_curve(y_test, y_pred_grd_lm)

# The gradient boosted model by itself
y_pred_grd = grd.predict_proba(X_test)[:, 1]
fpr_grd, tpr_grd, _ = roc_curve(y_test, y_pred_grd)

# The random forest model by itself
#y_pred_rf = rf.predict_proba(X_test)[:, 1]
#fpr_rf, tpr_rf, _ = roc_curve(y_test, y_pred_rf)

plt.figure(1)
plt.plot([0, 1], [0, 1], 'k--')
#plt.plot(fpr_rt_lm, tpr_rt_lm, label='RT + LR')
#plt.plot(fpr_rf, tpr_rf, label='RF')
#plt.plot(fpr_rf_lm, tpr_rf_lm, label='RF + LR')
plt.plot(fpr_grd, tpr_grd, label='GBT')
plt.plot(fpr_grd_lm, tpr_grd_lm, label='GBT + LR')
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.title('ROC curve')
plt.legend(loc='best')
plt.show()
