import useSWR from 'swr';
import NGL from './NGLComponent';
import PDB from './PDBCard';

import ExpansionPanel from '@material-ui/core/ExpansionPanel';
import ExpansionPanelSummary from '@material-ui/core/ExpansionPanelSummary';
import ExpansionPanelDetails from '@material-ui/core/ExpansionPanelDetails';
import ExpandMoreIcon from '@material-ui/icons/ExpandMore';

const fetcher = url => fetch(url).then(r => r.json())

const MoleculeComponent = ({ id, smiles }) => {
  const { data, error } = useSWR(`/api/smiles/${smiles}`, fetcher, { refreshInterval: 0 });
  if (error) return <div>failed to load</div>
  if (!data) return <div>loading...</div>

  let { chembl, pubchem } = data;
  console.log(data.bindings);

  return(
    <ExpansionPanel defaultExpanded>
      <ExpansionPanelSummary expandIcon={<ExpandMoreIcon />} id={ "ep" + id } >
        <p>SMILES: { smiles }</p>
      </ExpansionPanelSummary>
      <ExpansionPanelDetails>
        <PDB data={pubchem} />
        <p><img src="https://www.ebi.ac.uk/chembl/api/data/image/${ chembl.id }.svg" height="200px" /></p>

        <h2>Binding data</h2>
        <NGL
          data={{ filename: "https://files.rcsb.org/download/4hhb.pdb" }}
          viewportId={ "viewport-" + id }
        />
      </ExpansionPanelDetails>
    </ExpansionPanel>
  );
};

export default MoleculeComponent;
//{data.bindings.map((structure, index) => (
//  <p key={structure}>{structure}</p>
//  ))}

//<span title="Download SDF" style="float: right;"><a href="https://files.rcsb.org/ligands/view/{{ data.pdb.id }}_model.sdf"><i className="ti-export"></i></a></span>