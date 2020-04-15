const fetch = require('node-fetch');
const Pubchem = require('../models/pubchem');

// https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial
class PubchemFactory {
    static URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/<SMILES>/JSON';

    static getProperty(props, label, name) {
        let i;
        if (label == "label") {
            i = props.findIndex(k => k.urn.label == name);
        }
        else {
            i = props.findIndex(k => k.urn.name == name);
        }
        if (i != -1) {
            if (props[i].value.fval) {
                return props[i].value.fval;
            }
            if (props[i].value.ival) {
                return props[i].value.ival;
            }
            if (props[i].value.sval) {
                return props[i].value.sval;
            }
        }
    }

    static parse(json) {
        if (!json.PC_Compounds) return;
        //console.log(json.PC_Compounds[0].props);
        let props = json.PC_Compounds[0].props;
        let cid = json.PC_Compounds[0].id.id.cid;
        let smiles = this.getProperty(props, 'label', 'SMILES');
        let name = this.getProperty(props, 'label', 'IUPAC Name');
        let formula = this.getProperty(props, 'label', 'Molecular Formula');
        let mass = this.getProperty(props, 'label', 'Mass');
        let logp = this.getProperty(props, 'label', "Log P");
        let hydrogen_bond_acceptor = this.getProperty(props, "name", "Hydrogen Bond Acceptor");
        let hydrogen_bond_donor = this.getProperty(props, "name", "Hydrogen Bond Donor");
        let rotatable_bond = this.getProperty(props, "name", "Rotatable Bond");
        let polar_surface_area = this.getProperty(props, "name", "Polar Surface Area");
        let heavy_atoms = json.PC_Compounds[0].count.heavy_atom;
        let atom_chiral = json.PC_Compounds[0].count.atom_chiral;
        let tautomers = json.PC_Compounds[0].count.tautomers;
        let pubchem = new Pubchem(cid, smiles, name, formula, mass, logp,
                    heavy_atoms, atom_chiral, tautomers, hydrogen_bond_acceptor,
                    hydrogen_bond_donor, rotatable_bond, polar_surface_area);
        return(pubchem);
    }

    static async fetch(smiles) {
        let url = this.URL.replace('<SMILES>', encodeURIComponent(smiles));
        try {
            const res = await fetch(url);
            let data = await res.json();
            return(this.parse(data));
        } catch (err) {
            console.log(err);
        }
    }
}

module.exports = PubchemFactory;
