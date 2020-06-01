import parse from 'urlencoded-body-parser';
import NGL from '../components/molecule/NGLComponent';

import Grid from '@material-ui/core/Grid';
import Paper from '@material-ui/core/Paper';
import Layout from '../components/layout/Layout';
import MoleculeComponent from '../components/molecule/MoleculeComponent';

let molecules = [];

const Molecules = ({ smiles }) => {
  return(
    <Layout>
      <Grid container spacing={3}>
        <Grid item xs={6}>
          <Paper>
            <h1>Results</h1>
            <NGL
              data={{ filename: "https://files.rcsb.org/download/4hhb.pdb" }}
              viewportId="viewport-1"
            />
            {smiles.map((option, index) => (
              <MoleculeComponent key={option} smiles={option} />
            ))}
          </Paper>
        </Grid>
      </Grid>
    </Layout>
  );
};

export async function getServerSideProps(context) {
  const query = await parse(context.req);
  const { structure, structuresList } = query;
  if (structure) {
    molecules = [structure];
  }
  else if (structuresList) {
    molecules = [...new Set(structuresList.match(/[^\r\n]+/g).map((item) => item.trim()))];
  }
  return { props: { smiles: molecules } };
};

export default Molecules;