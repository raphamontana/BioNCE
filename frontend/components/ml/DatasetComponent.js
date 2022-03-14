import React from 'react';
import { Box, Button, Grid, TextField, Typography } from '@material-ui/core';
import { ExpansionPanel, ExpansionPanelDetails, ExpansionPanelSummary } from './ExpansionPanel';
import { Add as AddIcon } from '@material-ui/icons';
import AlertDialog from './AlertDialog';
import JSMEComponent from '../draw/JSMEComponent';
import StructuresUploadButton from '../draw/StructuresUploadButton';
import StructuresListComponent from '../draw/StructuresListComponent';
import useStyles from '../layout/style';
import SMILES from '../../libs/models/SMILES';

const DatasetComponent = ({ dataset, setDataset, handleNext }) => {
  const classes = useStyles();
  const [open, setOpen] = React.useState(false);
  const [expanded, setExpanded] = React.useState('panel1');
  const [structuresTextField, setStructuresTextField] = React.useState('');
  const [jsmeSmiles, setJSMESmiles] = React.useState('');

  const submit = () => {
    if (dataset.length !== 0) {
      handleNext();
    }
    else {
      setOpen(true);
    }
  }

  const addToDataset = (newStructures) => {
    let arr = newStructures.match(/[^\r\n]+/g);
    if (arr !== null) {
      arr = arr.map((item) => item.trim().split(' ')[0]);
      for (let i = 0; i < arr.length; i++) {
        if (!SMILES.isValid(arr[i])) { arr.splice(i, 1); }
      };
      const structuresSet = new Set(dataset);
      arr.forEach(e => structuresSet.add(e));
      if (structuresSet.size > dataset.length) {
        setDataset([...structuresSet]);
      }
    }
  };

  const handleExpansion = (panel) => (event, newExpanded) => {
    setExpanded(newExpanded ? panel : false);
  };

  return (
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
              <Typography>#1 Draw</Typography>
            </ExpansionPanelSummary>
            <ExpansionPanelDetails>
              <Grid
                container
                direction="column"
                justify="flex-start"
              >
                <Grid item>
                  <JSMEComponent callback={s => setJSMESmiles(s)} />
                  <Box display="inline" alignItems="center">
                    <TextField
                      id="structure"
                      name="structure"
                      label="SMILES"
                      variant="outlined"
                      className={classes.textfield}
                      value={jsmeSmiles}
                      onChange={e => setJSMESmiles(e.target.value)}
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
                </Grid>
              </Grid>
            </ExpansionPanelDetails>
          </ExpansionPanel>
          <ExpansionPanel square expanded={expanded === 'panel2'} onChange={handleExpansion('panel2')}>
            <ExpansionPanelSummary aria-controls="panel2d-content" id="panel2d-header">
              <Typography>#2 Write</Typography>
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
                    rows={15}
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
          <ExpansionPanel square expanded={expanded === 'panel3'} onChange={handleExpansion('panel3')}>
            <ExpansionPanelSummary aria-controls="panel3d-content" id="panel3d-header">
              <Typography>#3 Upload file</Typography>
            </ExpansionPanelSummary>
            <ExpansionPanelDetails>
              <StructuresUploadButton addStructures={addToDataset} />
            </ExpansionPanelDetails>
          </ExpansionPanel>
        </Grid>
        <Grid item xs>
          <StructuresListComponent structures={dataset} setStructures={setDataset} />
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