/**
 * PDB Controller module.
 * @module controller/pdb
 */
const fetch = require('node-fetch');
const xml2js = require('xml2js');
const PDB = require('../models/pdb');

// https://www.rcsb.org/pdb/software/rest.do#smiles
// https://www.rcsb.org/structure/4klb
// http://www.rcsb.org/pdb/rest/customReport.xml?pdbids=1OXR&customReportColumns=structureId,structureTitle,ligandId,Doi
// https://www.bindingdb.org/bind/tabLuceneResult.jsp?thisInput=cruzipain&submit=Go
// http://bindingmoad.org/pdbrecords/index/1OXR
// https://www.ebi.ac.uk/pdbe/entry/search/index/?searchParams=%7B%22q_pdb_id%22:%5B%7B%22value%22:%221oxr%22,%22condition1%22:%22AND%22,%22condition2%22:%22Contains%22%7D%5D,%22resultState%22:%7B%22tabIndex%22:0,%22paginationIndex%22:1,%22perPage%22:%2210%22,%22sortBy%22:%22Sort%20by%22%7D%7D
// https://www.ebi.ac.uk/pdbe/graph-api/pdbe_doc/
// https://www.npmjs.com/package/pdbmine

/**
 * @class Chemical Structure searches with SMILES strings.
 */
class PDBFactory {
  static SEARCH_URL = "https://www.rcsb.org/pdb/rest/smilesQuery?smiles=<SMILES>&search_type=exact";
  static DESCRIBE_URL = "http://www.rcsb.org/pdb/rest/describePDB?structureId=<STRUCTUREID>";


  /**
   * Cast a XML file to a PDBStructure object.
   * 
   * @param {string} xml - The data from the URL.
   * @returns {PDBStructure} The PDBStructure object.
   */
  static parseStructure(xml) {
    if (!xml) return;
    // Parse the string to XML Object.
    let json;
    xml2js.parseString(xml, (_err, result) => {
      json = result;
    });
    // Get structure data.
    let structureId = json.PDBdescription.PDB[0].$.structureId;
    let description = json.PDBdescription.PDB[0].$.title;
    let resolution = json.PDBdescription.PDB[0].$.resolution;
    // Create the object.
    let structure = new PDB.PDBStructure(structureId, description, resolution);
    return(structure);
  }


  /**
   * Return data given a structure ID.
   * 
   * @async
   * @param {string} structureId - The protein ID to fetch.
   * @returns {PDBStructure} A protein structure.
   */
  static async fetchStructure(structureId) {
    let url = this.DESCRIBE_URL.replace('<STRUCTUREID>', structureId);
    try {
      const res = await fetch(url);
      let xml = await res.text();
      let structure = this.parseStructure(xml);
      return(structure);
    } catch (err) {
      console.log(err);
    }
  }


  /**
   * Cast a XML file to a PDBLigand object.
   * 
   * @param {string} xml - The data from the URL.
   * @returns {PDBLigand} The PDBLigand object.
   */
  static async parseLigand(xml) {
    if (!xml) return;
    // Parse the string to XML Object.
    let json;
    xml2js.parseString(xml, (_err, result) => {
      json = result;
    });
    // Get ligand data.
    let chemicalID = json.smilesQueryResult.ligandInfo[0].ligand[0].$.chemicalID;
    let molecularWeight = json.smilesQueryResult.ligandInfo[0].ligand[0].$.molecularWeight;
    let chemicalName = json.smilesQueryResult.ligandInfo[0].ligand[0].chemicalName;
    let formula = json.smilesQueryResult.ligandInfo[0].ligand[0].formula;
    let smiles = json.smilesQueryResult.ligandInfo[0].ligand[0].smiles;
    // Get structures list.
    let structures = [];
    for (let l in json.smilesQueryResult.ligandInfo[0].ligand) {
      let structureId = json.smilesQueryResult.ligandInfo[0].ligand[l].$.structureId;
      let structure = await this.fetchStructure(structureId);
      structures.push(structure);
    }
    // Create the object.
    let ligand = new PDB.PDBLigand(chemicalID, chemicalName, formula, smiles, molecularWeight, structures);
    return(ligand);
  }


  /**
   * Search for ligands and PDB IDs based on a SMILES query.
   * 
   * @async
   * @param {string} smiles - The query.
   * @returns {Promise<Object>} A ligand object.
   */
  static async searchLigandBySmiles(smiles) {
    let url = this.SEARCH_URL.replace('<SMILES>', encodeURIComponent(smiles));
    try {
      const res = await fetch(url);
      if (!res.ok) return;
      let xml = await res.text();
      let pdb = await this.parseLigand(xml);
      //console.log(pdb);
      return(pdb);
    } catch (err) {
      console.log(err);
    }
  }
}


module.exports = PDBFactory;
