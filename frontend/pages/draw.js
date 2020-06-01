import { useState } from 'react';
import Box from '@material-ui/core/Box';
import Button from '@material-ui/core/Button';
import Grid from '@material-ui/core/Grid';
import Paper from '@material-ui/core/Paper';
import TextField from '@material-ui/core/TextField';

import AddIcon from '@material-ui/icons/Add';
import ClearIcon from '@material-ui/icons/Clear';
import SearchIcon from '@material-ui/icons/Search';
import UploadIcon from '@material-ui/icons/CloudUpload';

import Layout from '../components/layout/Layout';
import useStyles from '../components/layout/style';
import LoadDrugButton from '../components/draw/LoadDrugButton';
import JSMEComponent, { JSMEAcknowledgement } from '../components/draw/JSMEComponent';

const Draw = () => {
  const classes = useStyles();
  const [smiles, setSmiles] = useState('');
  const [smilesList, setSmilesList] = useState('');

  const concatSmilesList = (structure) => {
    if (structure !== '') {
      if (smilesList !== '') {
        setSmilesList(smilesList + '\n' + structure);
      }
      else {
        setSmilesList(structure);
      }
    }
  }

  const handleCopy = () => {
    if (smiles !== '') {
      if (smilesList !== '') setSmilesList(smilesList + '\n' + smiles);
      else setSmilesList(smiles);
    }
  }

  return(
    <Layout>
      <Grid container spacing={2}>
        <Grid item xs>
          <Paper className={`${classes.paper} ${classes.jsme}`} >
            Draw smiles
            <JSMEComponent callBack={e => setSmiles( e )} />
            <form action="/molecules" method="POST">
              <Box component="div" display="inline">
                <TextField
                  id="structure"
                  label="SMILES"
                  variant="outlined"
                  value={smiles}
                  onChange={e => setSmiles( e.target.value )}
                  className={classes.textfield}
                />
                  <Button
                    variant="contained"
                    color="primary"
                    startIcon={<AddIcon />}
                    onClick={handleCopy}
                    className={classes.button}
                  >
                    Copy
                  </Button>
                  <Button
                    type="submit"
                    variant="contained"
                    color="primary"
                    startIcon={<SearchIcon />}
                    className={classes.button}
                  >
                    Search
                  </Button>
              </Box>
            </form>
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
                type="submit"
                variant="contained"
                color="primary"
                startIcon={<UploadIcon />}
              >
                Upload
              </Button>
            </label>
          </Paper>
        </Grid>

        <Grid item xs>
          <Paper className={classes.paper}>
            <form action="/molecules" method="POST">
              <TextField
                id="structuresList"
                variant="outlined"
                label="Enter a list of SMILES here:"
                value={smilesList}
                onChange={e => setSmilesList(e.target.value)}
                multiline
                rows={18}
                className={classes.textfield}
              />
              <Box component="div" display="inline">
                <Button
                  variant="outlined"
                  color="secondary"
                  className={classes.button}
                  startIcon={<ClearIcon />}
                  onClick={() => setSmilesList('')}
                >
                  Clear
                </Button>
                <LoadDrugButton setSmilesList={concatSmilesList} />
                <Button
                  type="submit"
                  variant="contained"
                  color="primary"
                  className={classes.button}
                  startIcon={<SearchIcon />}
                >
                  Search
                </Button>
              </Box>
            </form>
          </Paper>
        </Grid>
      </Grid>
      <JSMEAcknowledgement />
    </Layout>
  );
}

export default Draw;