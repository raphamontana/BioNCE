import { makeStyles } from '@material-ui/styles';
import { Avatar, Button, Card, CardActions, CardContent, CardHeader, FormControl,
  IconButton, InputLabel, MenuItem , Select, Tooltip } from '@material-ui/core';

import NGL from './NGLComponent';

const useStyles = makeStyles(() => ({
  card: {
    minWidth: 300,
    maxWidth: 345,
  }
}));

const BindingCard = ({ data }) => {
  const classes = useStyles();

  return(
    <Card className={classes.card}>
      <NGL
        data={{ filename: "https://files.rcsb.org/download/4hhb.pdb" }}
        viewportId={ "viewport-" + data }
      />
    </Card>
  );
}

export default BindingCard;