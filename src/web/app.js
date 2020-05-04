const compression = require('compression');
const createError = require('http-errors');
const express = require('express');
const helmet = require('helmet')
const logger = require('morgan');
const nunjucks = require('nunjucks')

// Middleware
const app = express();
app.use(compression());
app.use(express.json()); // for parsing application/json
app.use(express.text());
app.use(express.urlencoded({
    extended: true
})) // for parsing application/x-www-form-urlencoded
app.use(express.static('public'));
app.use(helmet());
app.use(logger("combined"))

nunjucks.configure('views', {
    autoescape: true,
    express: app
});

// Routes
app.use('/', require('./routes/index'));
app.use('/dbsearch', require('./routes/dbsearch'));
app.use('/molecule', require('./routes/molecule'));
app.use('/smiles', require('./routes/smiles'));

// Catch 404 and forward to error handler.
app.use(function (req, res, next) {
    next(createError(404));
});

// Error handler.
app.use(function (err, req, res, next) {
    // Set locals, only providing error in development.
    res.locals.message = err.message;
    res.locals.error = req.app.get('env') === 'development' ? err : {};

    // Render the error page.
    res.status(err.status || 500);
    res.render('base.html');
});

let port = 80;
app.listen(port, () => console.log("'App listening on port http://localhost:" + port + ' .'));