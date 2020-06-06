import useSWR from 'swr';

import ExpansionPanel from '@material-ui/core/ExpansionPanel';
import ExpansionPanelDetails from '@material-ui/core/ExpansionPanelDetails';
import ExpansionPanelSummary from '@material-ui/core/ExpansionPanelSummary';
import ExpandMoreIcon from '@material-ui/icons/ExpandMore';

import BindingCard from './BindingCard';
import ChEMBLCard from './ChEMBLCard';
import PubChemCard from './PubChemCard';
import PDBCard from './PDBCard';

const fetcher = url => fetch(url).then(r => r.json())

const MoleculeComponent = ({ id, smiles }) => {
  const { data, error } = useSWR(`/api/smiles/${smiles}`, fetcher, { refreshInterval: 0 });
  if (error) return <div>failed to load</div>
  if (!data) return <div>loading...</div>

  let { chembl, pubchem, pdb } = data;
  console.log(data.bindings);

  return(
    <ExpansionPanel defaultExpanded>
      <ExpansionPanelSummary expandIcon={<ExpandMoreIcon />} id={ "ep" + id } >
        <p>SMILES: { smiles }</p>
      </ExpansionPanelSummary>
      <ExpansionPanelDetails>
        <PubChemCard data={pubchem} />
        <ChEMBLCard data={chembl} />
        <PDBCard data={pdb} />
        <BindingCard data={pdb.structures[0].ligands[0]} />
      </ExpansionPanelDetails>
    </ExpansionPanel>
  );
};

export default MoleculeComponent;
//<span title="Download SDF" style="float: right;"><a href="https://files.rcsb.org/ligands/view/{{ data.pdb.id }}_model.sdf"><i className="ti-export"></i></a></span>