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
            <Typography variant="h4">
              Contact Us
            </Typography>
            <Typography variant="body2">
              We would love to hear from you!<br />
              Have comments, suggestions or corrections?<br />
              Send us an email at <a href="mailto:bionce@nequimed.usp.br">bionce@nequimed.usp.br</a>
            </Typography>
            <Typography variant="h6">
              Feedback Form
            </Typography>
            <form>
              <Typography component="h2" variant="h6" color="primary" gutterBottom>
                Name
                First Name
                Last Name
                Email
                Found what I needed.
                Yes, I found what I needed.
                No, I did not find what I needed.
                What were you looking for?

                Optional Survey
                Yes No
                Is the website easy to use?
                Yes
                No
                Are the help functions clear and concise?
                Yes
                No
                Are the FAQ's sufficient?
                Yes
                No
                Did the search tool function as expected?
                Yes
                No
                Did the filter tool function as expected?
                Yes
                No
                Suggestions to help us make it better.
              </Typography>
              <button type="submit">Submit</button>
            </form>
          </Paper>
        </Grid>
      </Grid>
    </Layout>
  );
};

export default Home;