
class ChEMBL {
    constructor(chembl_id, smiles, molformula, black_box_warning, chirality,
                max_phase, alogp, aromatic_rings, cx_logd, cx_logp,
                cx_most_apka, full_mwt, hba, hbd, heavy_atoms, psa,
                qed_weighted, ro3_pass, rtb, pref_name) {
        this.chembl_id = chembl_id;
        this.smiles = smiles;
        this.molformula = molformula;
        this.black_box_warning = black_box_warning;
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

module.exports = ChEMBL;
