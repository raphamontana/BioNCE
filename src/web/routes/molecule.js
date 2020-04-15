const express = require('express');
const router = express.Router();
const PubchemController = require('../controllers/pubchemController');
const ChemblController = require('../controllers/chemblController');
const PDBController = require('../controllers/pdbController');

router.get('/', async function (req, res) {
  let smiles = req.query.smiles;
  let pubchem = PubchemController.fetch(smiles);
  let chembl = ChemblController.fetch(smiles);
  let pdb = PDBController.searchLigandBySmiles(smiles);
  res.render('molecule.html', {smiles: smiles, pubchem: await pubchem,
                               chembl: await chembl, error: null});
})

router.post('/', async function (req, res) {
  let smilesList = req.body.smiles;
  smilesList = smilesList.match(/[^\r\n]+/g).map(item => item.trim());
  console.log(smilesList);
  let molecules = [];
  await Promise.all(smilesList.map(async (smiles) => {
    let pubchem = PubchemController.fetch(smiles);
    let chembl = ChemblController.fetch(smiles);
    let pdb = PDBController.searchLigandBySmiles(smiles);
    mol = { "smiles": smiles, "pubchem": await pubchem, "chembl": await chembl, "pdb": await pdb };
    molecules.push(mol);
  }));
  res.render('newmolecule.html', {"mols": molecules, "error": null});
})

module.exports = router;
