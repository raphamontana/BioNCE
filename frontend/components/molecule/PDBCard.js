import Image from 'material-ui-image'

import { makeStyles } from '@material-ui/styles';
import Avatar from '@material-ui/core/Avatar';
import Button from '@material-ui/core/Button';
import Card from '@material-ui/core/Card';
import CardActions from '@material-ui/core/CardActions';
import CardContent from '@material-ui/core/CardContent';
import CardHeader from '@material-ui/core/CardHeader';
import Grid from '@material-ui/core/Grid';
import IconButton from '@material-ui/core/IconButton';
import Tooltip from '@material-ui/core/Tooltip';
import LinkIcon from '@material-ui/icons/Link';

import MolecularFormula from './MolecularFormula';

const useStyles = makeStyles(() => ({
  card: {
    maxWidth: 345,
  }
}));

const Attribute = ({ name, value }) => {
  return(
    <Grid container justify="space-between">
      <Grid item>
        <Tooltip title={ name }>
          <span>{ name }:</span>
        </Tooltip>
      </Grid>
      <Grid item>
        { value }
      </Grid>
    </Grid>
  );
}

const PDBCard = ({ data }) => {
  const classes = useStyles();
  return(
    <Card className={classes.card}>
      <CardHeader
        avatar={
          <Avatar variant="rounded" alt="PubChem icon" src="https://pubchem.ncbi.nlm.nih.gov/pcfe/favicon/favicon.ico" />
        }
        action={
          <IconButton aria-label="PubChem link"
            href={ "https://pubchem.ncbi.nlm.nih.gov/compound/" + data.cid }
            target="_blank"
          >
            <LinkIcon />
          </IconButton>
        }
        title={ `Compound CID: ${data.cid}` }
        subheader="PubChem"
      />
      <CardContent>
        <Image style={{ width:'50%', height:'50%' }} width='100px' height='100px' alt="Compound structure"  src={`https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid=${data.cid}`} />
        <p align="center">

        </p>
        <Attribute name="Name" value={ data.name } />
        <MolecularFormula formula={ data.formula } />
        <Attribute name="Molecular Weight" value={ `${parseFloat(data.mass).toFixed(2)}g/mol` } />
        <Attribute name="LogP" value={ data.logp } />
        <Attribute name="Heavy atoms" value={ data.heavy_atoms } />
        <Attribute name="Atom chiral" value={ data.atom_chiral } />
        <Attribute name="Tautomers" value={ data.tautomers } />
        <Attribute name="Hydrogen Bond Acceptor" value={ data.hydrogen_bond_acceptor } />
        <Attribute name="Hydrogen Bond Donor" value={ data.hydrogen_bond_donor } />
        <Attribute name="Rotatable Bond" value={ data.rotatable_bond } />
        <Attribute name="Polar surface area" value={ data.polar_surface_area } />
        <CardActions>
          <Tooltip title="Fingerprint Tanimoto-based 2-dimensional similarity search.">
            <Button size="small" color="primary">
              Similar Structures Search
            </Button>
          </Tooltip>
          <Tooltip title="Standard substructure search, finds structures in the database that contain the input structure as a part.">
            <Button size="small" color="primary">
              Substructure Search
            </Button>
          </Tooltip>
        </CardActions>
      </CardContent>
    </Card>
  );
};

export default PDBCard;