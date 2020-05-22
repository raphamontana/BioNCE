import parse from 'urlencoded-body-parser';

import Grid from '@material-ui/core/Grid';
import Paper from '@material-ui/core/Paper';
import Layout from '../components/Layout/Layout';

const Smiles = ({ smiles }) => {
  return(
    <Layout>
      <Grid container spacing={3}>
        <Grid item xs={6}>
          <Paper>
            <h1>SMILES: { smiles }</h1>
            <p>Information</p>
          </Paper>
        </Grid>
      </Grid>
    </Layout>
  );
};

export async function getServerSideProps(context) {
  const { smiles } = await parse(context.req);
  if (typeof smiles !== "undefined") {
    return { props: { smiles } };
  }
  else {
    return { props: { smiles: null } };
  }
};

export default Smiles;