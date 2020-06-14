import { Button, Paper, TextField, Typography } from '@material-ui/core';
import useStyles from '../layout/style';

const DatasetComponent = ({ handleNext }) => {
  const classes = useStyles();
  const [dataset, setDataset] = React.useState('');

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
        onClick={handleNext}
        className={classes.button}
      >
        Next
      </Button>
    </>
  );
};

export default DatasetComponent;