import MoleculeController from "../../../libs/controllers/moleculeController";

export default async (req, res) => {
  const { query: { smiles } } = req;

  let molecule = await MoleculeController.searchBySmiles(smiles);

  res.status(200).json(molecule);
}