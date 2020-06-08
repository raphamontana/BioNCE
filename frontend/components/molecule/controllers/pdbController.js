/**
 * PDB Controller module.
 * @module controller/pdb
 */
const fetch = require("node-fetch");
const xml2js = require("xml2js");
const PDB = require("../models/pdb");

// https://www.rcsb.org/pdb/software/rest.do#smiles
// https://www.rcsb.org/structure/4klb
// http://www.rcsb.org/pdb/rest/customReport.xml?pdbids=1OXR&customReportColumns=structureId,structureTitle,ligandId,Doi
// https://www.bindingdb.org/bind/tabLuceneResult.jsp?thisInput=cruzipain&submit=Go
// http://bindingmoad.org/pdbrecords/index/1OXR
// https://www.ebi.ac.uk/pdbe/entry/search/index/?searchParams=%7B%22q_pdb_id%22:%5B%7B%22value%22:%221oxr%22,%22condition1%22:%22AND%22,%22condition2%22:%22Contains%22%7D%5D,%22resultState%22:%7B%22tabIndex%22:0,%22paginationIndex%22:1,%22perPage%22:%2210%22,%22sortBy%22:%22Sort%20by%22%7D%7D
// https://www.ebi.ac.uk/pdbe/graph-api/pdbe_doc/
// https://www.npmjs.com/package/pdbmine


function compare(a, b) {
  const structureA = a.id.toUpperCase();
  const structureB = b.id.toUpperCase();

  let comparison = 0;
  if (structureA > structureB) {
    comparison = 1;
  } else if (structureA < structureB) {
    comparison = -1;
  }
  return comparison;
}


/**
 * @class Chemical Structure searches with SMILES strings.
 */
class PDBFactory {
  static SMILES_URL =
    "https://www.rcsb.org/pdb/rest/smilesQuery?smiles=<SMILES>&search_type=exact";
  static STRUCTURES_URL =
    "http://www.rcsb.org/pdb/rest/describePDB?structureId=<STRUCTUREID>";
  static LIGANDS_URL =
    "https://www.rcsb.org/pdb/rest/ligandInfo?structureId=<STRUCTUREID>";

  static get RESOLUTION_THRESHOLD() {
    return 2;
  }

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
    let structure = {
      structureId,
      description,
      resolution
    };
    return structure;
  }


  static parseLigandInfo(xml) {
    if (!xml) return;
    // Parse the string to XML Object.
    let json;
    xml2js.parseString(xml, (_err, result) => {
      json = result;
    });
    // Get structure data.
    let ligandList = [];
    json.structureId.ligandInfo[0].ligand.forEach(ligand => {
      ligandList.push(ligand.$.chemicalID);
    });
    // Return the object.
    return (ligandList);
  }


  /**
   * Return data given a structure ID.
   *
   * @async
   * @param {string} structureId - The protein ID to fetch.
   * @returns {PDBStructure} A protein structure.
   */
  static async fetchStructure(structureId) {
    let structureUrl = this.STRUCTURES_URL.replace("<STRUCTUREID>", structureId);
    let ligandUrl = this.LIGANDS_URL.replace("<STRUCTUREID>", structureId);
    try {
      const structureRes = await fetch(structureUrl);
      if (!structureRes.ok) return;
      let structureXml = await structureRes.text();
      let structure = this.parseStructure(structureXml);
      const ligandRes = await fetch(ligandUrl);
      if (!ligandRes.ok) return;
      let ligandXML = await ligandRes.text();
      let ligands = this.parseLigandInfo(ligandXML);
      let pdbStructure = new PDB.PDBStructure(structure.structureId, structure.description, structure.resolution, ligands);
      return (pdbStructure);
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
  static async parsePDB(xml) {
    if (!xml) return;
    // Parse the string to XML Object.
    let json;
    xml2js.parseString(xml, (_err, result) => {
      json = result;
    });
    // Get ligand data.
    let ligand = json.smilesQueryResult.ligandInfo[0].ligand;
    let pdbLigand;
    if (ligand) {
      let chemicalID = ligand[0].$.chemicalID;
      let molecularWeight = ligand[0].$.molecularWeight;
      let chemicalName = ligand[0].chemicalName;
      let formula = ligand[0].formula;
      let smiles = ligand[0].smiles;
      // Get structures list.
      let structures = [];
      await Promise.all(
        ligand.map(async (element) => {
          let structure = await this.fetchStructure(element.$.structureId);
          structures.push(structure);
        })
      );
      structures.sort(compare);
      pdbLigand = new PDB.PDBLigand(chemicalID, chemicalName, formula, smiles, molecularWeight, structures);
    }
    // Create the object.
    return (pdbLigand);
  }


  /**
   * Search for ligands and PDB IDs based on a SMILES query.
   *
   * @async
   * @param {string} smiles - The query.
   * @returns {Promise<Object>} A ligand object.
   */
  static async searchLigandBySmiles(smiles) {
    /// TODO: Treat caracter '@' conversion.
    let url = this.SMILES_URL.replace("<SMILES>", encodeURIComponent(smiles));
    try {
      const res = await fetch(url);
      if (!res.ok) return;
      let xml = await res.text();
      return (await this.parsePDB(xml));
    } catch (err) {
      console.log(err);
    }
  }
}

module.exports = PDBFactory;