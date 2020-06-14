import Grid from '@material-ui/core/Grid';
import Paper from '@material-ui/core/Paper';
import Typography from '@material-ui/core/Typography';
import Layout from '../components/layout/Layout';
import Link from '../components/layout/Link';
import useStyles from '../components/layout/style';

const Help = () => {
  const classes = useStyles();

  return(
    <Layout>
      <Grid container spacing={3}>
        <Grid item xs={12}>
          <Paper className={`${classes.paper} ${classes.center}`} >
            <p><img src="images/logo.svg" alt="logo" width="512" /></p>
            <Typography component="h2" variant="h6" color="primary" gutterBottom>
              From a New Chemical Entity to a Bioactive New Chemical Entity.
            </Typography>
            <Typography component="body2">
              Welcome to BioNCE, an Intelligent Digital System by <a href="http://nequimed.iqsc.usp.br/">NEQUIMED/IQSC/USP</a>
            </Typography>
            <Typography variant="h6">
              How Do I Get Started?
            </Typography>
            <Typography variant="h6">
              Searching
            </Typography>
            <Typography variant="h6">
              Data sources:
            </Typography>
            <Typography variant="body2">
              PubChem, PDB, ChEMBL, BindingMOAD
            </Typography>
            <Typography variant="h6">
              <Link href="/contact">Get in touch</Link>
            </Typography>
            <Typography variant="h6">
              Useful links:
            </Typography>
            <Typography variant="body2">
              <a href="http://partridgejiang.github.io/Kekule.js/">Kekule.js</a><br />
              <a href="http://aquaria.ws/">Aquaria.ws</a><br />
              <a href="http://molview.org/">MolView</a><br />
              <a href="http://avogadro.cc/">Avogadro</a><br />
              <a href="http://openbabel.org/">OpenBabel</a><br />
              <a href="https://www.molsoft.com/moledit.html">Moledit</a><br />
              <a href="http://www.sciencegeek.net/Chemistry/chemware/chemware.shtml">Chemware</a><br />
              <a href="http://www.openmolecules.org/index.html">OpenMolecules</a><br />
              <a href="http://www.cheminfo.org/Chemistry/Cheminformatics/OpenChemLib_js/index.html">Cheminfo</a><br />
              <a href="http://www.cheminfo.org/ML/Regression/index.html">Cheminfo: Regression</a><br />
              <a href="https://www.ucd.ie/chem/chemint/index.html">Chemically intelligent online tutorials</a><br />
              2D chemical structure image recognition<br />
              <a href="https://lifescience.opensource.epam.com/imago/index.html">Imago</a><br />
              <a href="https://cactus.nci.nih.gov/cgi-bin/osra/index.cgi">OSRA</a><br />
              pubchem Related compounds download 7. bioactivity<br />
              ChEMBL similarity and substructure download
            </Typography>
          </Paper>
        </Grid>
      </Grid>
    </Layout>
  );
};

export default Help;