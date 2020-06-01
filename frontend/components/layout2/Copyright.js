import Box from '@material-ui/core/Box';
import MuiLink from '@material-ui/core/Link';
import Typography from '@material-ui/core/Typography';

export default function Copyright() {
  return (
    <Box pt={4}>
      <Typography variant="body2" color="textSecondary" align="center">
        {'© ' + new Date().getFullYear() + ' '}
        <MuiLink href="http://nequimed.iqsc.usp.br/" target="_blank" color="inherit">
          NEQUIMED/IQSC/USP
        </MuiLink>
        {'. All rights reserved.'}
      </Typography>
    </Box>
  );
}