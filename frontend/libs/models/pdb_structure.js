/**
 * @class A Protein Data Bank structure.
 */
class PDBStructure {
    /**
     * @constructs
     * @param {string} structureId - The ID of the structure.
     * @param {string} description - The chemical description.
     * @param {number} resolution - The resolution.
     * @param {string[]} ligands - The ligands.
     */
    constructor(id, description, resolution, ligands) {
        this.id = id;
        this.description = description;
        this.resolution = resolution;
        this.ligands = ligands;
    }
}

export default PDBStructure;