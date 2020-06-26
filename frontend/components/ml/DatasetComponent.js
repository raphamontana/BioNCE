import { withStyles } from '@material-ui/core/styles';
import { Button, Grid, IconButton, TextField, Typography, Paper } from '@material-ui/core';
import MuiExpansionPanel from '@material-ui/core/ExpansionPanel';
import MuiExpansionPanelSummary from '@material-ui/core/ExpansionPanelSummary';
import MuiExpansionPanelDetails from '@material-ui/core/ExpansionPanelDetails';
import AddIcon from '@material-ui/icons/Add';

import AlertDialog from "./AlertDialog"
import StructuresUploadButton from "../draw/StructuresUploadButton"
import JSMEComponent from "../draw/JSMEComponent"
import StructuresListComponent from "../draw/StructuresListComponent"

import useStyles from '../layout/style';


const ExpansionPanel = withStyles({
  root: {
    border: '1px solid rgba(0, 0, 0, .125)',
    boxShadow: 'none',
    '&:not(:last-child)': {
      borderBottom: 0,
    },
    '&:before': {
      display: 'none',
    },
    '&$expanded': {
      margin: 'auto',
    },
  },
  expanded: {},
})(MuiExpansionPanel);

const ExpansionPanelSummary = withStyles({
  root: {
    backgroundColor: 'rgba(0, 0, 0, .03)',
    borderBottom: '1px solid rgba(0, 0, 0, .125)',
    marginBottom: -1,
    minHeight: 56,
    '&$expanded': {
      minHeight: 56,
    },
  },
  content: {
    '&$expanded': {
      margin: '12px 0',
    },
  },
  expanded: {},
})(MuiExpansionPanelSummary);

const ExpansionPanelDetails = withStyles((theme) => ({
  root: {
    padding: theme.spacing(2),
  },
}))(MuiExpansionPanelDetails);


const DatasetComponent = ({ dataset, setDataset, handleNext }) => {
  const classes = useStyles();
  const [open, setOpen] = React.useState(false);
  const [expanded, setExpanded] = React.useState('panel1');
  const [structuresTextField, setStructuresTextField] = React.useState('');

  const submit = () => {
    if (dataset.length !== 0) {
      handleNext();
    }
    else {
      setOpen(true);
    }
  }

  const addToDataset = (structures) => {
    let structuresList = structures.match(/[^\r\n]+/g);
    if (structuresList !== null) {
      structuresList = structuresList.map((item) => item.trim().split(' ')[0]);
      const structuresSet = new Set(dataset);
      structuresList.forEach(e => structuresSet.add(e));
      if (structuresSet.size > dataset.length) {
        setDataset([...structuresSet]);
      }
    }
  };

  const handleExpansion = (panel) => (event, newExpanded) => {
    setExpanded(newExpanded ? panel : false);
  };

  return(
    <>
      <Grid
        container
        direction="row"
        justify="flex-start"
        alignItems="stretch"
        spacing={2}
      >
        <Grid item xs>
          <Typography variant="h6" color="initial">DATASET INPUT</Typography>
          <ExpansionPanel square expanded={expanded === 'panel1'} onChange={handleExpansion('panel1')}>
            <ExpansionPanelSummary aria-controls="panel1d-content" id="panel1d-header">
              <Typography>#1 Write</Typography>
            </ExpansionPanelSummary>
            <ExpansionPanelDetails>
              <Grid
                container
                direction="column"
                justify="flex-start"
              >
                <Grid item>
                  <TextField
                    id="dataset"
                    label="One SMILES per line."
                    value={structuresTextField}
                    onChange={(e) => setStructuresTextField(e.target.value)}
                    variant="outlined"
                    multiline
                    rows={10}
                    className={classes.textfield}
                  />
                </Grid>
                <Grid item>
                  <Button
                    variant="contained"
                    color="primary"
                    startIcon={<AddIcon />}
                    onClick={() => addToDataset(structuresTextField)}
                  >
                    Add
                  </Button>
                </Grid>
              </Grid>
            </ExpansionPanelDetails>
          </ExpansionPanel>
          <ExpansionPanel square expanded={expanded === 'panel2'} onChange={handleExpansion('panel2')}>
            <ExpansionPanelSummary aria-controls="panel2d-content" id="panel2d-header">
              <Typography>#2 Draw</Typography>
            </ExpansionPanelSummary>
            <ExpansionPanelDetails>
              <Paper>
                <JSMEComponent />
              </Paper>
            </ExpansionPanelDetails>
          </ExpansionPanel>
          <ExpansionPanel square expanded={expanded === 'panel3'} onChange={handleExpansion('panel3')}>
            <ExpansionPanelSummary aria-controls="panel3d-content" id="panel3d-header">
              <Typography>#3 Upload file</Typography>
            </ExpansionPanelSummary>
            <ExpansionPanelDetails>
              <StructuresUploadButton addStructures={addToDataset}/>
            </ExpansionPanelDetails>
          </ExpansionPanel>
        </Grid>
        <Grid item xs>
          <StructuresListComponent structures={dataset} setStructures={setDataset}/>
        </Grid>
      </Grid>
      <Button disabled={true} className={classes.button}>
        Back
      </Button>
      <Button
        variant="contained"
        color="primary"
        onClick={submit}
        className={classes.button}
      >
        Next
      </Button>
      <AlertDialog
        title="No input error."
        content="Add at least one structure to the dataset to query the model."
        open={open}
        setOpen={setOpen}
      />
    </>
  );
};

export default DatasetComponent;