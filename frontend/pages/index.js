import { SimpleImg } from 'react-simple-img';
import Grid from '@material-ui/core/Grid';
import Paper from '@material-ui/core/Paper';
import Typography from '@material-ui/core/Typography';
import Layout from '../components/layout/Layout';
import useStyles from '../components/layout/style';

const Home = () => {
  const classes = useStyles();

  return(
    <Layout>
      <Grid container spacing={3}>
        <Grid item xs={12}>
          <Paper className={`${classes.paper} ${classes.center}`} >
            <SimpleImg
              src="images/logo.svg"
              alt="logo"
              height="114px"
            />
            <Typography component="h2" variant="h6" color="primary" gutterBottom>
              From a New Chemical Entity to a Bioactive New Chemical Entity.
            </Typography>

            <div className="grid">
              <a href="/draw" className="card">
                <h3>Search &rarr;</h3>
                <p>Draw a structure and search the database.</p>
              </a>
              <a href="/dataset" className="card">
                <h3>Create a dataset &rarr;</h3>
                <p>Select molecules and features to download file.</p>
              </a>

              <a href="/ml" className="card">
                <h3>Discover &rarr;</h3>
                <p>Machine Learning to discover Bioactive New Chemical Entities.</p>
              </a>

              <a href="/documentation" className="card">
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
          </Paper>
        </Grid>
      </Grid>
    </Layout>
  );
};

export default Home;