import { Button, Grid, List, Paper, Typography } from '@material-ui/core';
import ResultItem from './ResultItem';
import useStyles from '../layout/style';

const ResultsComponent = ({ dataset, model, handleReset }) => {
  const classes = useStyles();

  return(
    <>
      <Paper className={classes.paper} >
        <Typography variant="h6" color="initial">RESULTS PRESENTATION ({model})</Typography>
        <Grid container spacing={2}>
          <Grid item >
            <List dense>
              {dataset.map((smiles, index) =>
                <ResultItem key={index} model={model} smiles={smiles}/>
              )}
            </List>
          </Grid>
        </Grid>
      </Paper>
      <Button onClick={handleReset} className={classes.button}>
        Reset
      </Button>
    </>
  );
}

export default ResultsComponent;