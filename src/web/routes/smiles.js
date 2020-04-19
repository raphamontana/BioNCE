const express = require('express');
const router = express.Router();

router.get('/', (req, res, next) => res.render('smiles.html'));

module.exports = router;