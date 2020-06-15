import { Button, Paper, TextField, Typography } from '@material-ui/core';
import AlertDialog from "./AlertDialog"
import useStyles from '../layout/style';

const DatasetComponent = ({ dataset, setDataset, handleNext }) => {
  const classes = useStyles();
  const [open, setOpen] = React.useState(false);

  const submitDataset = () => {
    if (dataset !== "") {
      handleNext();
    }
    else {
      setOpen(true);
    }
  }

  return(
    <>
      <Paper className={classes.paper} >
        <Typography variant="h6" color="initial">DATASET INPUT</Typography>
        <TextField
        id="dataset"
        label=""
        value={dataset}
        onChange={(e) => setDataset(e.target.value)}
        variant="outlined"
        multiline
        rows={10}
        className={classes.textfield}
        />
      </Paper>
      <Button disabled={1} className={classes.button}>
        Back
      </Button>
      <Button
        variant="contained"
        color="primary"
        onClick={submitDataset}
        className={classes.button}
      >
        Next
      </Button>
      <AlertDialog
        title="No input written"
        content="Insert the dataset to query the model."
        open={open}
        setOpen={setOpen}
      />
    </>
  );
};

export default DatasetComponent;