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
          </Paper>
        </Grid>
      </Grid>
    </Layout>
  );
};

export default Home;