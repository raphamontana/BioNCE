export default (req, res) => {
  const {
    query: { smiles },
  } = req
  smiles = ${smiles};
  res.statusCode = 200
  res.setHeader('Content-Type', 'application/json')

  res.end(JSON.stringify({ mol: smiles }))
}