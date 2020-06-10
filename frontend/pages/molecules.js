import parse from 'urlencoded-body-parser';
import Grid from '@material-ui/core/Grid';
import Layout from '../components/layout/Layout';
import MoleculeComponent from '../components/molecule/MoleculeComponent';

let molecules = [];

const Molecules = ({ smiles }) => {
  return(
    <Layout>
      <Grid container spacing={3}>
        <Grid item xs>
          {smiles.map((option, index) => (
            <MoleculeComponent key={index} id={index} smiles={option} />
          ))}
        </Grid>
      </Grid>
    </Layout>
  );
};

export async function getServerSideProps(context) {
  const query = await parse(context.req);
  const { structure, structuresList, structuresFile } = query;
  if (structure) {
    molecules = [structure];
  }
  else if (structuresList) {
    molecules = [...new Set(structuresList.match(/[^\r\n]+/g).map((item) => item.trim()))];
  }
  else if (structuresFile) {
    molecules = [...new Set(structuresFile.split(','))];
  }
  return { props: { smiles: molecules } };
};

export default Molecules;