const express = require('express');
const router = express.Router();

router.route('/')
  .get(function(req, res) {
    res.render('dbsearch.html');
  })
  .post(function(req, res) {
    res.render('dbsearch.html');
  });

module.exports = router;
