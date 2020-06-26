import { Box, Button, Grid, Paper, TextField } from '@material-ui/core';
import { Add as AddIcon,
         Clear as ClearIcon,
         Search as SearchIcon } from '@material-ui/icons';
import Layout from '../components/layout/Layout';
import useStyles from '../components/layout/style';
import LoadDrugButton from '../components/draw/LoadDrugButton';
import JSMEComponent, { JSMEAcknowledgement } from '../components/draw/JSMEComponent';
import StructuresUploadButton from '../components/draw/StructuresUploadButton';
import StructuresListComponent from '../components/draw/StructuresListComponent';


const Draw = () => {
  const classes = useStyles();
  const [jsmeSmiles, setJSMESmiles] = React.useState('');
  const [smilesList, setSmilesList] = React.useState('');
  const [dataset, setDataset] = React.useState([]);

  const addToDataset = (newStructures) => {
    let newStructuresList = newStructures.match(/[^\r\n]+/g);
    if (newStructuresList !== null) {
      newStructuresList = newStructuresList.map((item) => item.trim().split(' ')[0]);
      const structuresSet = new Set(dataset);
      newStructuresList.forEach(e => structuresSet.add(e));
      if (structuresSet.size > dataset.length) {
        setDataset([...structuresSet]);
      }
    }
  };

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

  const handleSubmit = (event) => {
    if (dataset.length === 0) {
      event.preventDefault();
    }
  }

  return(
    <Layout>
      <Grid container spacing={2}>
        <Grid item xs>
          <Paper className={`${classes.paper} ${classes.jsme}`} >
            Draw smiles
            <JSMEComponent callback={e => setJSMESmiles( e )} />
            <Box display="inline" alignItems="center">
              <TextField
                id="structure"
                name="structure"
                label="SMILES"
                variant="outlined"
                className={classes.textfield}
                value={jsmeSmiles}
                onChange={e => setJSMESmiles( e.target.value )}
              />
              <Button
                variant="contained"
                color="primary"
                className={classes.button}
                startIcon={<AddIcon />}
                onClick={() => addToDataset(jsmeSmiles)}
              >
                Add
              </Button>
            </Box>
          </Paper>

          <Paper className={classes.paper}>
            <TextField
              id="structuresList"
              name="structuresList"
              variant="outlined"
              label="One SMILES per line."
              value={smilesList}
              onChange={e => setSmilesList(e.target.value)}
              multiline
              rows={10}
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
                variant="contained"
                color="primary"
                className={classes.button}
                startIcon={<AddIcon />}
                onClick={() => addToDataset(smilesList)}
              >
                Add
              </Button>
            </Box>
          </Paper>

          <Paper className={classes.paper} style={{align: 'center'}}>
            <StructuresUploadButton addStructures={addToDataset}/>
            <p>* File format: Files that contain multiple SMILES codes need to list one code per line followed by the name of the compound (separated by a space or tab). If no compound name is given, 'unnamed' will be assigned. These files have to carry a '.smi' or '.smiles' extension.</p>
          </Paper>
        </Grid>

        <Grid item xs>
          <StructuresListComponent structures={dataset} setStructures={setDataset} />
          <form action="/molecules" method="POST" onSubmit={handleSubmit}>
            <input type="hidden" id="dataset" name="dataset" value={dataset} />
            <Button
              type="submit"
              variant="contained"
              color="primary"
              startIcon={<SearchIcon />}
              className={classes.button}
            >
              Search
            </Button>
          </form>
        </Grid>
      </Grid>
      <JSMEAcknowledgement />
    </Layout>
  );
}

export default Draw;