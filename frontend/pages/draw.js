import { useState } from 'react';
import { Box, Button, Grid, Paper, TextField } from '@material-ui/core';
import { Add as AddIcon,
         Clear as ClearIcon,
         Search as SearchIcon,
         CloudUpload as UploadIcon } from '@material-ui/icons';
import Layout from '../components/layout/Layout';
import useStyles from '../components/layout/style';
import LoadDrugButton from '../components/draw/LoadDrugButton';
import JSMEComponent, { JSMEAcknowledgement } from '../components/draw/JSMEComponent';
import StructuresUploadButton from '../components/draw/StructuresUploadButton';

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
            <JSMEComponent callback={e => setSmiles( e )} />
            <form action="/molecules" method="POST">
              <Box display="inline" alignItems="center">
                <TextField
                  id="structure"
                  name="structure"
                  label="SMILES"
                  variant="outlined"
                  className={classes.textfield}
                  value={smiles}
                  onChange={e => setSmiles( e.target.value )}
                />
                <Button
                  variant="outlined"
                  color="primary"
                  className={classes.button}
                  startIcon={<AddIcon />}
                  onClick={handleCopy}
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
            <StructuresUploadButton />
          </Paper>
        </Grid>

        <Grid item xs>
          <Paper className={classes.paper}>
            <form action="/molecules" method="POST">
              <TextField
                id="structuresList"
                name="structuresList"
                variant="outlined"
                label="Enter a list of SMILES here:"
                value={smilesList}
                onChange={e => setSmilesList(e.target.value)}
                multiline
                rows={18}
                className={classes.textfield}
              />
              <Box display="inline">
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