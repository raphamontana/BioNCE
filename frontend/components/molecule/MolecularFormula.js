// TODO Substituir esta biblioteca para formatar a formula molecular.
import HTMLParser from 'fast-html-parser';

import Grid from '@mui/material/Grid';
import Tooltip from '@mui/material/Tooltip';

const MolecularFormula = ({ formula }) => {
  let formattedFormula = formula.split('').map( (s) => {
    return s.match(/[0-9]/i) ? `<sub>${s}</sub>` : s;
  }).reduce((str, l) => str + l);

  return(
    <Grid container justifyContent="space-between">
      <Grid item>
        <Tooltip title="Molecular Formula" enterDelay={500}>
          <b>Molecular Formula:</b>
        </Tooltip>
      </Grid>
      <Grid item>
        { formula }
      </Grid>
    </Grid>
  );
}

export default MolecularFormula;
