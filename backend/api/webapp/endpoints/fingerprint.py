from rdkit import Chem
from rdkit.Chem import AllChem


def getFingerprint(smi):
    mol = Chem.MolFromSmiles(smi)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    return fp1
