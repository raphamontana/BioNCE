// Protein Data Bank

/**
 * @class A protein structure.
 */
class PDBStructure {
    /**
     * @constructs
     * @param {string} structureId - The ID of the ligand.
     * @param {string} description - The chemical name.
     * @param {number} resolution - The molecular formula.
     */
    constructor(structureId, description, resolution) {
        this.structureId = structureId;
        this.description = description;
        this.resolution = resolution;
        //this.doi = doi;
    }
}

/**
 * @class A protein ligand.
 */
class PDBLigand {
    /**
     * @constructs
     * @param {string} chemicalID - The ID of the ligand.
     * @param {string} chemicalName - The chemical name.
     * @param {string} formula - The molecular formula.
     * @param {string} smiles - The SMILES representation.
     * @param {number} molecularWeight - The molecular weight.
     * @param {Object[]} structures - A list of structures.
     */
    constructor(chemicalID, chemicalName, formula, smiles, molecularWeight, structures) {
        this.chemicalID = chemicalID;
        this.name = chemicalName;
        this.formula = formula;
        this.smiles = smiles;
        this.molecularWeight = molecularWeight;
        this.structures = structures;
    }
}

module.exports = {
  PDBLigand : PDBLigand,
  PDBStructure : PDBStructure
};
