import { Button, FormControl, InputLabel, MenuItem,
         Paper, Select, Typography } from '@material-ui/core';
import AlertDialog from "./AlertDialog"
import useStyles from '../layout/style';

const ModelSelectComponent = ({ model, setModel, handleBack, handleNext }) => {
  const classes = useStyles();
  const [open, setOpen] = React.useState(false);

  const models = {
    "covid-19-rf-v20.6": "COVID-19 (RF v20.6)",
    //"covid-19-fda-v20.6": "COVID-19 (FDA v20.6)",
    //"covid-19-rfe-v20.6": "COVID-19 (RFE v20.6)",
    "covid-19-mlp-v20.6": "COVID-19 (MLP v20.6)",
    "covid-19-gpc-v20.6": "COVID-19 (GPC v20.6)",
    "catl-catb-gpc-v20.6": "CatL-CatB (GPC v20.6)",
    "catl-catb-mlp-v20.6": "CatL-CatB (MLP v20.6)",
    "catl-catb-rf-v20.6": "CatL-CatB (RF v20.6)",
    //"catl-catb-fda-v20.6": "CatL-CatB (FDA v20.6)",
    //"catl-catb-rfe-v20.6": "CatL-CatB (RFE v20.6)",
    "catl-cats-gpc-v20.6": "CatL-CatS (GPC v20.6)",
    "catl-cats-mlp-v20.6": "CatL-CatS (MLP v20.6)",
    "catl-cats-rf-v20.6": "CatL-CatS (RF v20.6)"
    //"catl-cats-fda-v20.6": "CatL-CatS (FDA v20.6)",
    //"catl-cats-rfe-v20.6": "CatL-CatS (RFE v20.6)",
    //"drugbank-rf-v20.6": "DrugBank (RF v20.6)"
  };

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
            {Object.keys(models).map((model, index) =>
              <MenuItem key={index} value={model}>{models[model]}</MenuItem>
            )}
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