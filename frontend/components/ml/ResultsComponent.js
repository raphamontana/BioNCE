import { Button, Paper, Typography } from '@material-ui/core';
import useStyles from '../layout/style';

const ResultsComponent = ({ handleReset }) => {
  const classes = useStyles();

  return(
    <>
      <Paper className={classes.paper} >
        <Typography variant="h6" color="initial">RESULTS PRESENTATION</Typography>
      </Paper>
      <Button onClick={handleReset} className={classes.button}>
        Reset
      </Button>
    </>
  );
}

export default ResultsComponent;