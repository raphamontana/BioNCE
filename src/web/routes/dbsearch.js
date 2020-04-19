const express = require('express');
const router = express.Router();

router.get('/', (req, res, next) => res.render('dbsearch.html'));

module.exports = router;