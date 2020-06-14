import fetch from "node-fetch";
import ChEMBL from "../models/chembl";

//https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services
class ChEMBLFactory {
  static URL_ID =
    "/chembl/api/data/molecule/<ID>";
  static URL_SEARCH =
    "https://www.ebi.ac.uk/chembl/api/data/molecule/search?q=<QUERY>";
  static URL_SIMILARITY =
    "https://www.ebi.ac.uk/chembl/api/data/similarity/<SMILES>/100.json";
  static URL_SUBSTRUCTURE =
    "https://www.ebi.ac.uk/chembl/api/data/substructure/<SMILES>.json";

  static parse(json) {
    if (!json) return;
    let i = json.molecules.findIndex(
      (m) =>
        m.molecule_hierarchy.molecule_chembl_id ==
        m.molecule_hierarchy.parent_chembl_id
    );
    let mol = json.molecules[i];
    let chembl_id = mol.molecule_chembl_id;
    let smiles = mol.molecule_structures.canonical_smiles;
    let molformula = mol.molecule_properties.full_molformula;
    let black_box_warning = mol.black_box_warning;
    let chirality = mol.chirality;
    let max_phase = mol.max_phase;
    let alogp = mol.molecule_properties.alogp;
    let aromatic_rings = mol.molecule_properties.aromatic_rings;
    let cx_logd = mol.molecule_properties.cx_logd;
    let cx_logp = mol.molecule_properties.cx_logp;
    let cx_most_apka = mol.molecule_properties.cx_most_apka;
    let full_mwt = mol.molecule_properties.full_mwt;
    let hba = mol.molecule_properties.hba;
    let hbd = mol.molecule_properties.hbd;
    let heavy_atoms = mol.molecule_properties.heavy_atoms;
    let psa = mol.molecule_properties.psa;
    let qed_weighted = mol.molecule_properties.qed_weighted;
    let ro3_pass = mol.molecule_properties.ro3_pass;
    let rtb = mol.molecule_properties.rtb;
    let pref_name = mol.pref_name;
    let chembl = new ChEMBL(
      chembl_id,
      smiles,
      molformula,
      black_box_warning,
      chirality,
      max_phase,
      alogp,
      aromatic_rings,
      cx_logd,
      cx_logp,
      cx_most_apka,
      full_mwt,
      hba,
      hbd,
      heavy_atoms,
      psa,
      qed_weighted,
      ro3_pass,
      rtb,
      pref_name
    );
    return chembl;
  }

  static async similarity(smiles) {
    const url = this.URL_SIMILARITY.replace("<SMILES>", encodeURIComponent(smiles));
    try {
      const res = await fetch(url);
      const data = await res.json();
      return this.parse(data);
    } catch (err) {
      console.log(err);
    }
  }

  static async getByID(id) {
    const url = this.URL_ID.replace("<ID>", encodeURIComponent(id));
  }

  static async searchByName(query) {
    const url = this.URL_SEARCH.replace("<QUERY>", encodeURIComponent(query));
  }
}

export default ChEMBLFactory;