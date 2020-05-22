export default (req, res) => {
  const {
    query: { rid },
  } = req
  const requestId = 1;
  res.end( requestId );
}