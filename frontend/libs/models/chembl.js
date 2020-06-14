/**
 * @class The ChEMBL class.
 */
class ChEMBL {
    /**
     * @constructs
     * @param {string} id - The ChEMBL ID of the structure.
     * @param {string} smiles - Canonical smiles, generated using pipeline pilot.
     * @param {string} formula - Molecular formula for the full compound (including any salt).
     * @param {number} bbw - Indicates that the drug has a black box warning.
     * @param {string} chirality - Shows whether a drug is dosed as a racemic mixture (0), single stereoisomer (1) or is an achiral molecule (2).
     * @param {string} max_phase - Maximum phase of development reached for the compound (4 = approved). Null where max phase has not yet been assigned.
     * @param {string} alogp - Calculated ALogP.
     * @param {number} aromatic_rings - Number of aromatic rings.
     * @param {number} cx_logd - The cx_logd.
     * @param {number} cx_logp - The cx_logp.
     * @param {number} cx_most_apka - The cx_most_apka.
     * @param {number} full_mwt - The full mwt.
     * @param {number} hba - The hba.
     * @param {number} hbd - The hbd.
     * @param {number} heavy_atoms - Number of heavy atoms.
     * @param {number} psa - The Polar surface area.
     * @param {string} qed_weighted - The qed weighted.
     * @param {string} ro3_pass - Indicates whether the compound passes the rule-of-three (mw < 300, logP < 3 etc).
     * @param {string} rtb - Number rotatable bonds.
     * @param {string} pref_name - Preferred name for the molecule
     */
    constructor(id, smiles, formula, bbw, chirality,
        max_phase, alogp, aromatic_rings, cx_logd, cx_logp, cx_most_apka,
        full_mwt, hba, hbd, heavy_atoms, psa, qed_weighted, ro3_pass, rtb,
        pref_name) {
        this.id = id;
        this.smiles = smiles;
        this.formula = formula;
        this.bbw = bbw;
        this.chirality = chirality;
        this.max_phase = max_phase;
        this.alogp = alogp;
        this.aromatic_rings = aromatic_rings;
        this.cx_logd = cx_logd;
        this.cx_logp = cx_logp;
        this.cx_most_apka = cx_most_apka;
        this.full_mwt = full_mwt;
        this.hba = hba;
        this.hbd = hbd;
        this.heavy_atoms = heavy_atoms;
        this.psa = psa;
        this.qed_weighted = qed_weighted;
        this.ro3_pass = ro3_pass;
        this.rtb = rtb;
        this.pref_name = pref_name;
    }
}

export default ChEMBL;