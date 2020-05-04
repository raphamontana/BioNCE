const express = require("express");
const router = express.Router();
const MoleculeController = require("../controllers/moleculeController");

// @TODO Replace server-side await with client-side fetch.
// https://www.xul.fr/en/html5/fetch.php
// https://stackoverflow.com/questions/49813883/nunjucks-async-rendering-with-keystone-view
// https://mozilla.github.io/nunjucks/api.html#asynchronous
router.get("/", async function (req, res) {
  let smiles = req.query.smiles;
  if (!smiles) {
    res.render("molecule.html", {
      mols: null,
      error: "No SMILES to search."
    });
  } else {
    let molecules = await MoleculeController.searchBySmiles(smiles);
    res.render("molecule.html", {
      mols: molecules,
      error: null
    });
  }
});

router.post("/", async function (req, res) {
  let smiles = req.body.smiles;
  if (!smiles) {
    res.render("molecule.html", {
      mols: null,
      error: "No SMILES to search."
    });
  } else {
    let molecules = await MoleculeController.searchBySmiles(smiles);
    res.render("molecule.html", {
      mols: molecules,
      error: null
    });
  }
});

router.get("/:smiles", async function (req, res) {
  let smiles = req.params.smiles;
  res.send('user ' + req.params.smiles);
});

module.exports = router;