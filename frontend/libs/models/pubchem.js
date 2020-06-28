/**
 * @class The PubChem class.
 */
class Pubchem {
  constructor(cid, smiles, name, formula, mass, logp, heavy_atoms, atom_chiral,
    tautomers, hydrogen_bond_acceptor, hydrogen_bond_donor,
    rotatable_bond, polar_surface_area) {
    this.cid = cid;
    this.smiles = smiles;
    this.name = name;
    this.formula = formula;
    this.mass = mass;
    this.logp = logp;
    this.heavy_atoms = heavy_atoms;
    this.atom_chiral = atom_chiral;
    this.tautomers = tautomers;
    this.hydrogen_bond_acceptor = hydrogen_bond_acceptor;
    this.hydrogen_bond_donor = hydrogen_bond_donor;
    this.rotatable_bond = rotatable_bond;
    this.polar_surface_area = polar_surface_area;
  }
}

export default Pubchem;