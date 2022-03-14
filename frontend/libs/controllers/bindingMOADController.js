import axios from 'axios';
//import fetch from "node-fetch";
import Papa from "papaparse";

function compare(a, b) {
  const chainA = a.chain.toUpperCase();
  const chainB = b.chain.toUpperCase();

  let comparison = 0;
  if (chainA > chainB) {
    comparison = 1;
  } else if (chainA < chainB) {
    comparison = -1;
  }
  return comparison;
}


// https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial
class BindingMOADController {
  static URL = "http://bindingmoad.org/files/csv/<STRUCTURE>.csv";

  static parse(data) {
    if (!data) return;
    try {
      let csv = Papa.parse(data);
      let bindings = [];
      let proteinID;
      csv.data.forEach(row => {
        if (row[2]) {
          proteinID = row[2];
        }
        if (row[5]) {
          let ligand = row[3].split(':')[0];
          let chain = row[3].split(':')[1];
          let residue = row[3].split(':')[2];
          let affinityMeasure = row[5];
          switch (affinityMeasure) {
            case 'ic50':
              affinityMeasure = "IC<sub>50</sub>";
              break;
            case 'Kd':
              affinityMeasure = "K<sub>d</sub>";
              break;
            case 'Ki':
              affinityMeasure = "K<sub>i</sub>";
              break;
            case 'Ka':
              affinityMeasure = "K<sub>a</sub>";
              break;
            default:
              break;
          }
          let relationType = row[6];
          let affinityValue = row[7];
          let affinityUnit = row[8];
          let smiles = row[9];
          bindings.push(
            { proteinID, ligand, chain, residue, affinityMeasure, relationType, affinityValue, affinityUnit, smiles }
          );
        }
      });
      return (bindings.sort(compare));
    } catch (err) {
      console.log(err);
    }
    return (data);
  }

  static async fetch(pdbId) {
    let url = this.URL.replace('<STRUCTURE>', pdbId.toLowerCase());
    try {
      const resp = await axios.get(url);
      //const res = await fetch(url);
      return (this.parse(resp.data));
    } catch (err) {
      console.log(err);
    }
  }
}

export default BindingMOADController;