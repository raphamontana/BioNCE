import useSWR from "swr";

import LinearProgress from '@mui/material/LinearProgress';
import Accordion from '@mui/material/Accordion';
import AccordionSummary from '@mui/material/AccordionSummary';
import AccordionDetails from '@mui/material/AccordionDetails';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';

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
    <Accordion defaultExpanded>
      <AccordionSummary expandIcon={<ExpandMoreIcon />} id={ "ep" + id } >
        <p>SMILES: { smiles }</p>
      </AccordionSummary>
      <AccordionDetails>
        <PubChemCard data={pubchem} />
        <ChEMBLCard data={chembl} />
        <PDBCard data={pdb} />
        <BindingCard ligand={pdb} bindings={bindings} />
      </AccordionDetails>
    </Accordion>
  );
};

export default MoleculeComponent;