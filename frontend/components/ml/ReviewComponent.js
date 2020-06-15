import { Button, Paper, Typography } from '@material-ui/core';
import useStyles from '../layout/style';

const ReviewComponent = ({ dataset, model, handleBack, handleNext }) => {
  const classes = useStyles();

  return(
    <>
      <Paper className={classes.paper} >
        <Typography variant="h6" color="initial">REVIEW</Typography>
        <Typography variant="body2" color="initial">Dataset: {dataset}</Typography>
        <Typography variant="body2" color="initial">Model: {model}</Typography>
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
        Run
      </Button>
    </>
  );
}

export default ReviewComponent;