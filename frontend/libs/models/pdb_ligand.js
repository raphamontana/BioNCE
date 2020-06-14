import PDBStructure from './pdb_structure';

/**
 * @class A Protein Data Bank ligand.
 */
class PDBLigand {
    /**
     * @constructs
     * @param {string} id - The chemical ID of the ligand.
     * @param {string} name - The chemical name.
     * @param {string} formula - The molecular formula.
     * @param {string} smiles - The SMILES representation.
     * @param {number} weight - The molecular weight.
     * @param {Object[]} structures - A list of structures.
     */
    constructor(id, name, formula, smiles, weight, structures) {
        this.id = id;
        this.name = name;
        this.formula = formula;
        this.smiles = smiles;
        this.weight = weight;
        this.structures = structures;
    }

    // Getter
    get structuresList() {
        let sList = [];
        this.structures.forEach(s => {
            sList.push(s.id);
        });
        return (sList);
    }
}

export default PDBLigand;