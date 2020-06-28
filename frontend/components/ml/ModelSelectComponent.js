import { Button, FormControl, InputLabel, MenuItem,
         Paper, Select, Typography } from '@material-ui/core';
import AlertDialog from "./AlertDialog"
import useStyles from '../layout/style';

const ModelSelectComponent = ({ model, setModel, handleBack, handleNext }) => {
  const classes = useStyles();
  const [open, setOpen] = React.useState(false);

  const hasSelected = () => {
    if (model !== "") {
      handleNext();
    }
    else {
      setOpen(true);
    }
  }

  return(
    <>
      <Paper className={classes.paper} >
        <Typography variant="h6" color="initial">MODEL SELECTION</Typography>
        <FormControl variant="outlined" className={classes.formControl} >
          <InputLabel id="model-select-label">Model</InputLabel>
          <Select
            labelId="model-select-label"
            id="model-select"
            value={model}
            onChange={(e) => setModel(e.target.value)}
            label="Model"
          >
            <MenuItem value={"COVID-19-RF-v20.6"}>COVID-19 (RF v20.6)</MenuItem>
            <MenuItem value={"COVID-19-GPC-v20.6"}>COVID-19 (GPC v20.6)</MenuItem>
            <MenuItem value={"Chagas-ANN-v20.6"}>Chagas disease (ANN v20.6)</MenuItem>
            <MenuItem value={"Leishmaniasis-GPC-v20.6"}>Leishmaniasis (GPC v20.6)</MenuItem>
          </Select>
        </FormControl>
        <Typography variant="h6" color="initial">* About this model</Typography>
      </Paper>
      <Button onClick={handleBack} className={classes.button}>
        Back
      </Button>
      <Button
        variant="contained"
        color="primary"
        onClick={hasSelected}
        className={classes.button}
      >
        Next
      </Button>
      <AlertDialog
        title="No model selected"
        content="Choose a ML model to proceed."
        open={open}
        setOpen={setOpen}
      />
    </>
  );
}

export default ModelSelectComponent;