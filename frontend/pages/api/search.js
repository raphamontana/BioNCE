export default (req, res) => {
  const {
    query: { query },
  } = req
  const requestId = 1;
  res.end( requestId );
}
//https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/aspirin/json
//https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/assay/p68/json?limit=8
//https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/gene/egfr/jsonp?limit=5
//https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/taxonomy/mouse/json?limit=5