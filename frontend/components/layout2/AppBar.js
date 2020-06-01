import AppBar from '@material-ui/core/AppBar';
import Avatar from '@material-ui/core/Avatar';
import IconButton from '@material-ui/core/IconButton';
import Link from '@material-ui/core/Link';
import Toolbar from '@material-ui/core/Toolbar';
import useStyles from './style';

export default function appbar() {
  const classes = useStyles();

  return (
    <AppBar position="absolute" className={classes.appBar}>
      <Toolbar paddingright="24">
        <Link href="/">
          <IconButton className={classes.menuButton} >
            <Avatar
              src={"/images/logomark.svg"}
              alt={"Bionce"}
            />
          </IconButton>
        </Link>
      </Toolbar>
    </AppBar>
  );
}