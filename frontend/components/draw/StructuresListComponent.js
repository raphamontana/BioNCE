import { makeStyles } from '@material-ui/core/styles';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemIcon from '@material-ui/core/ListItemIcon';
import ListItemSecondaryAction from '@material-ui/core/ListItemSecondaryAction';
import ListItemText from '@material-ui/core/ListItemText';
import Checkbox from '@material-ui/core/Checkbox';
import IconButton from '@material-ui/core/IconButton';
import Paper from '@material-ui/core/Paper';
import Typography from '@material-ui/core/Typography';
import DeleteIcon from '@material-ui/icons/Delete';

const useStyles = makeStyles((theme) => ({
  list: {
    width: '100%',
    maxWidth: '100%',
    backgroundColor: theme.palette.background.paper,
    position: 'relative',
    overflow: 'auto',
    height: 300,
    maxHeight: 300,
  },
  paper: {
    display: 'flex',
    flexDirection: 'column',
    padding: theme.spacing(2),
  },
}));


const StructuresListComponent = ({ structures, setStructures }) => {
  const classes = useStyles();
  const [checked, setChecked] = React.useState([]);
  //const [smiles, setSmiles] = React.useState(["CC(=O)OC1=CC=CC=C1C(=O)O", "CCCCCC", "H2O", "SO4"]);

  const handleToggle = (value) => () => {
    const currentIndex = checked.indexOf(value);
    const newChecked = [...checked];
    if (currentIndex === -1) {
      newChecked.push(value);
    } else {
      newChecked.splice(currentIndex, 1);
    }
    setChecked(newChecked);
  };

  const handleDelete = (value) => () => {
    let pos = checked.indexOf(value);
    if (pos !== -1) {
      const newChecked = [...checked];
      newChecked.splice(pos, 1);
      setChecked(newChecked);
    }
    pos = structures.indexOf(value);
    const newStructures = [...structures];
    newStructures.splice(pos, 1);
    setStructures(newStructures);
  };

  return (
    <Paper className={classes.paper}>
      <Typography variant="h6" color="initial">STRUCTURES LIST</Typography>
      <List className={classes.list}>
        {structures.length === 0 ?
          <Typography variant="body1" color="initial">Empty</Typography>
          : null
        }
        {structures.map((value, index) => {
          const labelId = `checkbox-list-label-${index}`;

          return(
            <ListItem key={index} role={undefined} dense button onClick={handleToggle(value)}>
              <ListItemIcon>
                <Checkbox
                  edge="start"
                  checked={checked.indexOf(value) !== -1}
                  tabIndex={-1}
                  disableRipple
                  inputProps={{ 'aria-labelledby': labelId }}
                />
              </ListItemIcon>
              <ListItemText id={labelId} primary={value} />
              <ListItemSecondaryAction>
                <IconButton edge="end" aria-label="comments" onClick={handleDelete(value)}>
                  <DeleteIcon />
                </IconButton>
              </ListItemSecondaryAction>
            </ListItem>
          );
        })}
      </List>
    </Paper>
  );
}

export default StructuresListComponent;