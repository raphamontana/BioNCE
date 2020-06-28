import { Molecule } from "openchemlib";

/**
 * @class The SMILES class.
 */
class SMILES {
  /**
   * @constructs
   * @param {string} id - The ID of the structure.
   */
  constructor(id) {
    this.id = id;
  }

  static isValid(smiles) {
    try {
      const molfile = Molecule.fromSmiles(smiles, {noCoordinates: true, noStereo: true});
      if (smiles === '' || molfile.toIsomericSmiles() === '') {
        return false;
      }
      else return true;
    }
    catch {
      return false;
    }
  }
}

export default SMILES;