import Grid from '@material-ui/core/Grid';
import Tooltip from '@material-ui/core/Tooltip';

const Attribute = ({ name, value, tooltip }) => {
  let tooltipText = tooltip ? tooltip : name;
  return(
    <Grid container justify="space-between">
      <Grid item>
        <Tooltip title={ tooltipText } enterDelay={ 500 }>
          <b>{ name }:</b>
        </Tooltip>
      </Grid>
      <Grid item>
        { value }
      </Grid>
    </Grid>
  );
}

export default Attribute;