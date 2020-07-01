import fetch from "node-fetch";

export default async (req, res) => {
  const { query: { slug } } = req;
  const model = slug[0]
  const smiles = slug[1];

  let active;
  let url = `https://bionce.herokuapp.com/?model=${encodeURIComponent(model)}&smiles=${encodeURIComponent(smiles)}`; //encodeURIComponent

  try {
    const res = await fetch(url);
    const data = await res.json();
    console.log(data);
    if (data.prediction === "'SMILES Parse Error'") {

    }
    if (data.prediction === '1') {
      if (model.includes("covid")) {
        active = "Active";
      }
      else active = "CatL";
    }
    else {
      if (model.includes("covid")) {
        active = "Inactive";
      }
      else active = "CatS";
    }
  } catch (err) {
    console.log(err);
    res.status(400).json({ err: err });
  }
  res.status(200).json({ active: active });
}