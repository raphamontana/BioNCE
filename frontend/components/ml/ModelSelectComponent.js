import { Button, FormControl, InputLabel, MenuItem,
         Paper, Select, Typography } from '@material-ui/core';
import useStyles from '../layout/style';

const ModelSelectComponent = ({ handleBack, handleNext }) => {
  const classes = useStyles();
  const [model, setModel] = React.useState('');

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
          <MenuItem value={"ANN"}>Artificial Neural Network</MenuItem>
          <MenuItem value={"GPC"}>Gaussian Process Classifier</MenuItem>
          <MenuItem value={"RF"}>Random Forests</MenuItem>
        </Select>
      </FormControl>
      </Paper>
      <Button onClick={handleBack} className={classes.button}>
        Back
      </Button>
      <Button
        variant="contained"
        color="primary"
        onClick={handleNext}
        className={classes.button}
      >
        Next
      </Button>
    </>
  );
}

export default ModelSelectComponent;