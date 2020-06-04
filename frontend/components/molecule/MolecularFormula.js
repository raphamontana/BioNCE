import ReactHtmlParser from 'react-html-parser';

import Grid from '@material-ui/core/Grid';
import Tooltip from '@material-ui/core/Tooltip';

const MolecularFormula = ({ formula }) => {
  let formattedFormula = formula.split('').map( (s) => {
    return s.match(/[0-9]/i) ? `<sub>${s}</sub>` : s;
  }).reduce((str, l) => str + l);
  return(
    <Grid container justify="space-between">
      <Grid item>
        <Tooltip title="Molecular Formula">
          <span>Molecular Formula:</span>
        </Tooltip>
      </Grid>
      <Grid item>
        { ReactHtmlParser(formattedFormula) }
      </Grid>
    </Grid>
  );
}

export default MolecularFormula;