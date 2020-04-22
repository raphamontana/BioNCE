const bodyParser = require('body-parser');
const compression = require('compression');
const createError = require('http-errors');
const express = require('express');
const helmet = require('helmet')
//const morgan = require("morgan")
const nunjucks = require('nunjucks')

const app = express();
app.use(bodyParser.urlencoded({
    extended: false
}));
app.use(bodyParser.json());
app.use(bodyParser.text());
app.use(compression());
app.use(express.static('public'));
app.use(helmet());
//app.use(morgan("combined"))

nunjucks.configure('views', {
    autoescape: true,
    express: app
});

// Views
app.use('/', require('./routes/index'));
app.use('/dbsearch', require('./routes/dbsearch'));
app.use('/molecule', require('./routes/molecule'));
app.use('/smiles', require('./routes/smiles'));

let port = 3002
app.listen(port, () => console.log("'App listening on port " + port + '.'));