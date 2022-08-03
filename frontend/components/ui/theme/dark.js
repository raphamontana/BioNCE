import { createTheme } from '@material-ui/core/styles'

const darkPalette = {
  palette: {
    type: 'dark',
  },
};
const themeName = 'BioNCE dark theme';

// Create a theme instance.
export default createTheme({ darkPalette, themeName });
