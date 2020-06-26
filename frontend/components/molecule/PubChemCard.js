import { SimpleImg } from 'react-simple-img';
import { makeStyles } from '@material-ui/styles';
import Avatar from '@material-ui/core/Avatar';
import Button from '@material-ui/core/Button';
import Card from '@material-ui/core/Card';
import CardActions from '@material-ui/core/CardActions';
import CardContent from '@material-ui/core/CardContent';
import CardHeader from '@material-ui/core/CardHeader';
import IconButton from '@material-ui/core/IconButton';
import Tooltip from '@material-ui/core/Tooltip';
import GetAppIcon from '@material-ui/icons/GetApp';

import MolecularFormula from './MolecularFormula';
import Attribute from './Attribute';

const useStyles = makeStyles((theme) => ({
  card: {
    minWidth: 290,
    maxWidth: 345,
    margin: theme.spacing(0.5),
  }
}));

const PubChemCard = ({ data }) => {
  const classes = useStyles();

  if (typeof data === "undefined") {
    return(null);
  }

  return(
    <Card className={classes.card}>
      <CardHeader
        avatar={
          <a
            aria-label="PubChem link"
            href={ `https://pubchem.ncbi.nlm.nih.gov/compound/${data.cid}` }
            target="_blank"
          >
            <Avatar variant="rounded" alt="PubChem icon" src="https://pubchem.ncbi.nlm.nih.gov/pcfe/favicon/favicon.ico" />
          </a>
        }
        action={
          <IconButton
            aria-label="PubChem link"
            href={`https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/${data.cid}/JSON/?response_type=save&response_basename=compound_CID_${data.cid}`}
          >
            <GetAppIcon />
          </IconButton>
        }
        title={ `Compound ${data.cid}` }
        subheader="PubChem"
      />
      <CardContent>
        <div align="center">
          <SimpleImg
            src={`https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid=${data.cid}`}
            alt="A two-dimensional representation of the compound"
            height={100}
          />
        </div>
        <Attribute name="Name" value={ data.name } />
        <MolecularFormula formula={ data.formula } />
        <Attribute name="Molecular Mass" value={ `${parseFloat(data.mass).toFixed(2)} g/mol` } />
        <Attribute name="LogP" value={ data.logp } />
        <Attribute name="Heavy atoms" value={ data.heavy_atoms } tooltip="Number of non-hydrogen atoms." />
        <Attribute name="Atom chiral" value={ data.atom_chiral } />
        <Attribute name="Tautomers" value={ data.tautomers } />
        <Attribute name="H-Bond Acceptor" value={ data.hydrogen_bond_acceptor } tooltip="Number of hydrogen-bond acceptors in the structure." />
        <Attribute name="H-Bond Donor" value={ data.hydrogen_bond_donor } tooltip="Number of hydrogen-bond donors in the structure." />
        <Attribute name="Rotatable Bond" value={ data.rotatable_bond } tooltip="Number of rotatable bonds." />
        <Attribute name="Polar surface area" value={ data.polar_surface_area } tooltip="Topological polar surface area, computed by the algorithm described in the paper by Ertl et al." />
        <CardActions>
          <Tooltip title="Fingerprint Tanimoto-based 2-dimensional similarity search." enterDelay={500}>
            <Button size="small" color="primary" href={`https://pubchem.ncbi.nlm.nih.gov/#query=CID${data.cid}%20structure&tab=similarity`}>
              Similar structures
            </Button>
          </Tooltip>
          <Tooltip title="Standard substructure search, finds structures in the database that contain the input structure as a part." enterDelay={500}>
            <Button size="small" color="primary">
              Substructure search
            </Button>
          </Tooltip>
        </CardActions>
      </CardContent>
    </Card>
  );
};

export default PubChemCard;