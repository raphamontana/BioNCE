import Router from 'next/router';
import PropTypes from 'prop-types';
import { withStyles } from '@material-ui/core/styles';
import Box from '@material-ui/core/Box';
import Button from '@material-ui/core/Button';
import ButtonGroup from '@material-ui/core/ButtonGroup';
import Grid from '@material-ui/core/Grid';
import Link from '@material-ui/core/Link';
import Paper from '@material-ui/core/Paper';
import TextField from '@material-ui/core/TextField';
import Typography from '@material-ui/core/Typography';

import AddIcon from '@material-ui/icons/Add';
import ClearIcon from '@material-ui/icons/Clear';
import SearchIcon from '@material-ui/icons/Search';
import UploadIcon from '@material-ui/icons/CloudUpload';

import LoadExampleButton from '../components/LoadExampleButton';
import Layout from '../components/Layout/Layout';
import JSMEComponent from '../components/JSMEComponent';

const useStyles = theme => ({
  paper: {
    padding: theme.spacing(2),
    display: 'flex',
    overflow: 'auto',
    flexDirection: 'column',
    minWidth: '369px',
    marginBottom: theme.spacing(3),
  },
  textfield: {
    display: 'flex',
    marginTop: theme.spacing(2),
  },
  textarea: {
    flexDirection: 'column',
    height: '436px',
  },
  button: {
    margin: theme.spacing(1),
  },
  input: {
    display: 'none',
  },
});


class Draw extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      smiles: "",
      smilesList: "",
    };
  }

  handleSubmit = () => {
    Router.push({
      pathname: '/smiles',
      query: { s: this.state.smiles },
    })
  }

  handleCopy = (e) => {
    if (this.state.smiles !== '') {
      if (this.state.smilesList !== '') {
        this.setState({smilesList: this.state.smilesList + '\n' + this.state.smiles});
      }
      else {
        this.setState({smilesList: this.state.smiles});
      }
    }
  }

  handleClear = (e) => {
    this.setState({ smilesList: '' });
  }

  handleLoadExample = (drug) => {
    let smiles;
    if (drug==="asa")  // Aspirin
      smiles = "CC(=O)OC1=CC=CC=C1C(=O)O";
    else if (drug==="ibu")  // Ibuprofen
      smiles = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O";
    else if (drug==="ome")  // Omeprazole
      smiles = "CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=CC(=C3)OC";
    else if (drug==="sil")  // Sildenafil
      smiles = "CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12";
    else if (drug==="vas")  // Atorvastatin
      smiles = "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4";
    else if (drug==="azi")  // Azithromycin
      smiles = "CC[C@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@H](N(C)C)[C@H]2O)[C@](C)(O)C[C@@H](C)CN(C)[C@H](C)[C@@H](O)[C@]1(C)O";
  }

  handleSmiles = (e) => {
    let mol = null;
    if (typeof e === "string") {
      mol = e;
    }
    else {
      mol = e.target.value;
    }
    this.setState({ smiles: mol });
    console.log(mol);
  }

  handleSmilesList = (e) => {
    let mols = e.target.value;
    this.setState({ smilesList: mols });
    console.log(mols);
  }

  render() {
    const { classes } = this.props;
    return(
      <>
        <Layout>
          <Grid container spacing={3}>
            <Grid item xs={6}>
              <Paper className={classes.paper}>
                Draw smiles
                <JSMEComponent callBack={this.handleSmiles} />
                <Box component="div" display="inline">
                  <TextField
                    label="SMILES"
                    variant="outlined"
                    value={this.state.smiles}
                    onChange={this.handleSmiles}
                    className={classes.textfield}
                  />
                  <Button
                    variant="contained"
                    color="primary"
                    startIcon={<SearchIcon />}
                    className={classes.button}
                    onClick={this.handleSubmit}
                  >
                    Search
                  </Button>
                  <Button
                    variant="contained"
                    color="primary"
                    startIcon={<AddIcon />}
                    className={classes.button}
                    onClick={this.handleCopy}
                  >
                    Copy
                  </Button>
                  <form action="/smiles" method="post">
                    <label htmlFor="fname">First name:</label>
                    <input type="text" id="fname" name="fname" /><br /><br />
                    <label htmlFor="lname">Last name:</label>
                    <input type="text" id="lname" name="lname" /><br /><br />
                    <input type="submit" value="Submit" />
                  </form>
                </Box>
              </Paper>

              <Paper className={classes.paper} style={{align: 'center'}}>
                <input
                  accept=".csv"
                  className={classes.input}
                  id="contained-button-file"
                  multiple
                  type="file"
                />
                <label htmlFor="contained-button-file">
                  <Button
                    variant="contained"
                    color="primary"
                    startIcon={<UploadIcon />}
                    onClick={this.handleSubmit}
                  >
                    Upload
                  </Button>
                </label>
              </Paper>
            </Grid>

            <Grid item xs={6}>
              <Paper className={classes.paper}>
                <TextField
                  variant="outlined"
                  label="Enter a list of SMILES here:"
                  value={this.state.smilesList}
                  onChange={this.handleSmilesList}
                  multiline
                  rows={18}
                />
                <ButtonGroup color="secondary" aria-label="outlined secondary button group">
                  <Button
                    variant="contained"
                    color="primary"
                    startIcon={<ClearIcon />}
                    className={classes.button}
                    onClick={this.handleClear}
                  >
                    Clear
                  </Button>
                  <LoadExampleButton />
                  <Button
                    variant="contained"
                    color="primary"
                    startIcon={<SearchIcon />}
                    className={classes.button}
                    onClick={this.handleSubmit}
                  >
                    Search
                  </Button>
                </ButtonGroup>
              </Paper>
            </Grid>
          </Grid>

          <Box pt={4}>
            <Typography variant="body2" color="textSecondary" align="center">
              <Link href="http://peter-ertl.com/jsme/" target="_blank">JSME</Link>
              : Bruno Bienfait and Peter Ertl,{' '}
              <Link href="http://www.jcheminf.com/content/5/1/24" target="_blank">
                JSME: a free molecule editor in JavaScript
              </Link>
              , <i>Journal of Cheminformatics</i> <strong>5</strong>:24 (2013).
            </Typography>
          </Box>
        </Layout>
      </>
    );
  }
}

Draw.propTypes = {
  classes: PropTypes.object.isRequired,
};

export default withStyles(useStyles)(Draw);