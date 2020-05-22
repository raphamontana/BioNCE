//import { makeStyles } from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';
import Paper from '@material-ui/core/Paper';
import Typography from '@material-ui/core/Typography';
import Layout from '../components/Layout/Layout';
import useStyles from '../components/Layout/style';

/*
const useStyles = makeStyles((theme) => ({
  paper: {
    padding: theme.spacing(2),
    display: 'flex',
    overflow: 'auto',
    flexDirection: 'column',
    alignItems: 'center',
  },
}));
*/

const Home = () => {
  const classes = useStyles();

  return(
    <Layout>
      <Grid container spacing={3}>
        <Grid item xs={12}>
          <Paper className={classes.paper} alignitems={"center"} >
            <p>Welcome to BioNCE, an Intelligent Digital System by <a href="http://nequimed.iqsc.usp.br/">NEQUIMED/IQSC/USP</a></p>
            <p><img src="images/logo.svg" alt="logo" width="512" /></p>
            <Typography component="h2" variant="h6" color="primary" gutterBottom>
              From a New Chemical Entity to a Bioactive New Chemical Entity.
            </Typography>

            <div className="grid">
              <a href="/draw" className="card">
                <h3>Search &rarr;</h3>
                <p>Draw a structure and search the database.</p>
              </a>
              <a href="/molecule" className="card">
                <h3>Create a database &rarr;</h3>
                <p>Select molecules and features to download file.</p>
              </a>

              <a href="/ml" className="card">
                <h3>Run Machine Learning &rarr;</h3>
                <p>Discover Bioactive New Chemical Entities using Machine Learning.</p>
              </a>

              <a href="/examples" className="card">
                <h3>Documentation &rarr;</h3>
                <p>Find information about BioNCE features.</p>
              </a>
            </div>
            <style jsx>{`
            .grid {
              display: flex;
              align-items: center;
              justify-content: center;
              flex-wrap: wrap;
              max-width: 800px;
              margin-top: 3rem;
            }
            .card {
              margin: 1rem;
              flex-basis: 45%;
              padding: 1.5rem;
              text-align: left;
              color: inherit;
              text-decoration: none;
              border: 1px solid #eaeaea;
              border-radius: 10px;
              transition: color 0.15s ease, border-color 0.15s ease;
            }
            .card:hover,
            .card:focus,
            .card:active {
              color: #0070f3;
              border-color: #0070f3;
            }
            .card h3 {
              margin: 0 0 1rem 0;
              font-size: 1.5rem;
            }
            .card p {
              margin: 0;
              font-size: 1.25rem;
              line-height: 1.5;
            }
            @media (max-width: 600px) {
              .grid {
                width: 100%;
                flex-direction: column;
              }
            }
            `}</style>

            <h4>Links interessantes:</h4>
            <p>
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
              <a href="https://cactus.nci.nih.gov/cgi-bin/osra/index.cgi">OSRA</a>
            </p>
          </Paper>
        </Grid>
      </Grid>
    </Layout>
  );
};

export default Home;