import useSWR from "swr";

import LinearProgress from "@material-ui/core/LinearProgress";
import ExpansionPanel from "@material-ui/core/ExpansionPanel";
import ExpansionPanelDetails from "@material-ui/core/ExpansionPanelDetails";
import ExpansionPanelSummary from "@material-ui/core/ExpansionPanelSummary";
import ExpandMoreIcon from "@material-ui/icons/ExpandMore";

import SMILES from "../../libs/models/SMILES";
import BindingCard from "./BindingCard";
import ChEMBLCard from "./ChEMBLCard";
import PubChemCard from "./PubChemCard";
import PDBCard from "./PDBCard";

const fetcher = url => fetch(url).then(r => r.json());

const MoleculeComponent = ({ id, smiles }) => {
  if (!SMILES.isValid(smiles)) {
    console.log("Not done yet");
  }
  const { data, error } = useSWR(`/api/smiles/${smiles}`, fetcher, { refreshInterval: 0 });
  if (error) return <div>failed to load</div>
  if (!data) return <LinearProgress />

  let { chembl, pubchem, pdb, bindings } = data;

  return(
    <ExpansionPanel defaultExpanded>
      <ExpansionPanelSummary expandIcon={<ExpandMoreIcon />} id={ "ep" + id } >
        <p>SMILES: { smiles }</p>
      </ExpansionPanelSummary>
      <ExpansionPanelDetails>
        <PubChemCard data={pubchem} />
        <ChEMBLCard data={chembl} />
        <PDBCard data={pdb} />
        <BindingCard ligand={pdb} bindings={bindings} />
      </ExpansionPanelDetails>
    </ExpansionPanel>
  );
};

export default MoleculeComponent;