import { useRouter } from 'next/router';
import Grid from '@material-ui/core/Grid';
import Paper from '@material-ui/core/Paper';
import Layout from '../components/Layout/Layout';
import { useEffect } from 'react';

const Smiles = ( {myvar} ) => {
  const router = useRouter()
  const { smiles } = router.query

  useEffect(() => {
    //window.history.pushState(null, null, "/smiles");
  });

  return (
    <Layout>
      <Grid container spacing={3}>
        <Grid item xs={6}>
          <Paper>
            <h1>Title</h1>
            <p>SMILES: {myvar}</p>
          </Paper>
        </Grid>
      </Grid>
    </Layout>
  )
}

export async function getServerSideProps(context) {
  const { query } = context.req;

  return {
    props: React.Component.getInitialProps ? await React.Component.getInitialProps(ctx): {
        query
    }
  }
}

export default Smiles;