import clsx from 'clsx';
import Avatar from '@material-ui/core/Avatar';
import Divider from '@material-ui/core/Divider';
import Drawer from '@material-ui/core/Drawer';
import HomeIcon from '@material-ui/icons/Home';
import IconButton from '@material-ui/core/IconButton';
import LibraryBooksIcon from '@material-ui/icons/LibraryBooks';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemAvatar from '@material-ui/core/ListItemAvatar';
import ListItemIcon from '@material-ui/core/ListItemIcon';
import ListItemText from '@material-ui/core/ListItemText';
import ListSubheader from '@material-ui/core/ListSubheader';
import LocalHospitalIcon from '@material-ui/icons/LocalHospital';
import MenuIcon from '@material-ui/icons/Menu';
import MuiLink from '@material-ui/core/Link';
import SearchIcon from '@material-ui/icons/Search';
import TimelineIcon from '@material-ui/icons/Timeline';
import NextLink from '../Link';
import useStyles from './style';

export default function appbar() {
  const classes = useStyles();

  const [open, setOpen] = React.useState(false);
  const handleDrawer = () => {
    if (open) {
      setOpen(false);
    }
    else setOpen(true);
  };

  return (
    <Drawer
        variant="permanent"
        classes={{
          paper: clsx(classes.drawerPaper, !open && classes.drawerPaperClose),
        }}
        open={open}
    >
      <div className={classes.toolbarIcon} />
      <Divider />
      <List>
      <NextLink href="/">
        <ListItem button>
          <ListItemIcon>
            <IconButton edge="start"
              color="inherit"
              aria-label="open drawer"
              onClick={handleDrawer}
              className={clsx(classes.menuButton)}>
              <MenuIcon />
            </IconButton>
          </ListItemIcon>
        </ListItem>
      </NextLink>
      <NextLink href="/">
        <ListItem button>
          <ListItemIcon>
            <HomeIcon />
          </ListItemIcon>
          <ListItemText primary="Home" />
        </ListItem>
      </NextLink>
      <NextLink href="/draw">
        <ListItem button>
          <ListItemIcon>
            <SearchIcon />
          </ListItemIcon>
          <ListItemText primary="Search" />
        </ListItem>
      </NextLink>
      <ListItem button>
        <ListItemIcon>
          <TimelineIcon />
        </ListItemIcon>
        <ListItemText primary="Molecule" />
      </ListItem>
      <ListItem button>
        <ListItemIcon>
          <LocalHospitalIcon />
        </ListItemIcon>
        <ListItemText primary="Diceases" />
      </ListItem>
      <ListItem button>
        <ListItemIcon>
          <LibraryBooksIcon />
        </ListItemIcon>
        <ListItemText primary="Documentation" />
      </ListItem>
      </List>
      <Divider />
      <List>
      <ListSubheader inset>Brought to you by:</ListSubheader>
        <MuiLink color="inherit" href="http://nequimed.iqsc.usp.br/" target="_blank">
          <ListItem button>
            <ListItemAvatar>
              <Avatar
                alt={"Nequimed"}
                src={"/images/nequimed-avatar.png"}
              />
            </ListItemAvatar>
            <ListItemText primary="Nequimed" />
          </ListItem>
        </MuiLink>
        <MuiLink color="inherit" href="http://www.crob.eesc.usp.br/" target="_blank">
          <ListItem button>
            <ListItemAvatar>
              <Avatar
                alt={"Centro de Robótica da USP"}
                src={"/images/crob-avatar.png"}
              />
            </ListItemAvatar>
            <ListItemText primary="CROB" />
          </ListItem>
        </MuiLink>
        <MuiLink color="inherit" href="http://www5.iqsc.usp.br/" target="_blank">
          <ListItem button>
            <ListItemAvatar>
              <Avatar
                alt={"Instituto de Química de São Carlos"}
                src={"/images/iqsc-usp-logo.png"}
              />
            </ListItemAvatar>
            <ListItemText primary="IQSC" />
          </ListItem>
        </MuiLink>
        <MuiLink color="inherit" href="https://prp.usp.br/" target="_blank">
          <ListItem button>
            <ListItemAvatar>
              <Avatar
                alt={"Pró-Reitoria de Pesquisa da USP"}
                src={"/images/prp-usp-avatar.png"}
              />
            </ListItemAvatar>
            <ListItemText primary="PRP" />
          </ListItem>
        </MuiLink>
        <MuiLink color="inherit" href="http://www.fapesp.br/" target="_blank">
          <ListItem button>
            <ListItemAvatar>
              <Avatar
                alt={"FAPESP"}
                src={"/images/fapesp-logo.png"}
              />
            </ListItemAvatar>
            <ListItemText primary="FAPESP" />
          </ListItem>
        </MuiLink>
      </List>
    </Drawer>
  );
}