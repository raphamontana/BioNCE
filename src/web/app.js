const bodyParser = require('body-parser')
const compression = require('compression');
const createError = require('http-errors');
const express = require('express');
const helmet = require('helmet')
const nunjucks = require('nunjucks')

const app = express()
app.use(bodyParser.urlencoded({ extended: false }))
app.use(bodyParser.json())
app.use(bodyParser.text())
app.use(compression());
app.use(express.static('public'));
app.use(helmet())

nunjucks.configure('views', {
    autoescape: true,
    express: app
});

// Views
app.use('/', require('./routes/index'));
app.use('/dbsearch', require('./routes/dbsearch'));
app.use('/molecule', require('./routes/molecule'));
app.use('/smiles', require('./routes/smiles'));

app.listen(80, () => console.log('App listening on port 80.'));
