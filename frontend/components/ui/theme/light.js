// material-ui/examples/nextjs/src/theme.js
import { createMuiTheme, responsiveFontSizes } from '@material-ui/core/styles';
import { red } from '@material-ui/core/colors';

const palette = {
  primary: { main: '#556cd6' },
  secondary: { main: '#19857b' },
  error: { main: red.A400 },
  background: { default: '#fff' },
};
const themeName = 'BioNCE light theme';

// Create a theme instance.
export default createMuiTheme({ palette, themeName });