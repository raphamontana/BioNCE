import BMOADController from "./bindingMOADController";
import ChemblController from "./chemblController";
import PDBController from "./pdbController";
import PubchemController from "./pubchemController";

class MoleculeController {

  static async searchBySmiles(smiles) {

  //const res = await fetch('https://.../posts');
  //const posts = await res.json()

    const pubchem = PubchemController.fetch(smiles);
    const chembl = ChemblController.fetch(smiles);
    const pdb = await PDBController.searchLigandBySmiles(smiles);
    const bindings = [];
    if (pdb) {
      await Promise.all(
        pdb.structuresList.map(async structure => {
          let binding = await BMOADController.fetch(structure);
          if (binding.length > 0) {
            bindings.push(binding);
          }
        })
      );
    }
    let molecule = {
      smiles: smiles,
      pubchem: await pubchem,
      chembl: await chembl,
      pdb: pdb,
      bindings: bindings
    };
    return(molecule);
  }

}

export default MoleculeController;