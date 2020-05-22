export default (req, res) => {
  const {
    query: { smiles },
  } = req;
  res.statusCode = 200;
  res.setHeader('Content-Type', 'application/json');
  res.end(JSON.stringify({ mol: smiles }));
}