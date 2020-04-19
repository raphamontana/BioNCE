const PubchemController = require("./pubchemController");
const ChemblController = require("./chemblController");
const PDBController = require("./pdbController");
const BMOADController = require("./bindingMOADController");

class MoleculeController {

    static async searchBySmiles(smiles) {
        let smilesList = smiles.match(/[^\r\n]+/g).map((item) => item.trim());
        let molecules = [];
        await Promise.all(
            smilesList.map(async smiles => {
                let pubchem = PubchemController.fetch(smiles);
                let chembl = ChemblController.fetch(smiles);
                let pdb = await PDBController.searchLigandBySmiles(smiles);
                let bindings = [];
                await Promise.all(
                    pdb.structuresList.map(async structure => {
                        let binding = await BMOADController.fetch(structure);
                        if (binding.length > 0) {
                            bindings.push(binding);
                        }
                    }));
                molecules.push({
                    smiles: smiles,
                    pubchem: await pubchem,
                    chembl: await chembl,
                    pdb: pdb,
                    bindings: bindings
                });
            })
        );
        return (molecules);
    }

}

module.exports = MoleculeController;