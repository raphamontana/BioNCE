import { SimpleImg } from 'react-simple-img';
import { makeStyles } from '@material-ui/styles';
import { Avatar, Card, CardContent, CardHeader, IconButton, Tooltip } from '@material-ui/core';
import GetAppIcon from '@material-ui/icons/GetApp';
import Attribute from './Attribute';
import MolecularFormula from './MolecularFormula';

const useStyles = makeStyles((theme) => ({
  card: {
    minWidth: 290,
    maxWidth: 345,
    margin: theme.spacing(0.5),
  }
}));

const ChEMBLCard = ({ data }) => {
  const classes = useStyles();

  if (typeof data === "undefined") {
    return (null);
  }

  return (
    <Card className={classes.card}>
      <CardHeader
        avatar={
          <a href={`https://www.ebi.ac.uk/chembl/compound_report_card/${data.id}/`} target="_blank">
            <Avatar variant="rounded" alt="ChEMBL icon" src="https://www.ebi.ac.uk/chembl/k8s/static/chembl/img/icons/chembl_logo_pink.png" />
          </a>
        }
        action={
          <Tooltip title="Download model file." enterDelay={500}>
            <IconButton
              aria-label="ChEMBL link"
              href={`https://www.ebi.ac.uk/chembl/api/data/molecule/${data.id}.sdf`}
              target="_blank"
            >
              <GetAppIcon />
            </IconButton>
          </Tooltip>
        }
        title={data.id}
        subheader="ChEMBL"
      />
      <CardContent>
        <div align="center">
          <SimpleImg
            src={`https://www.ebi.ac.uk/chembl/api/data/image/${data.id}.svg`}
            alt="Compound structure"
            height={100}
          />
        </div>
        <Attribute name="Name" value={data.pref_name} />
        <MolecularFormula formula={data.formula} />
        <Attribute name="Molecular Mass" value={`${parseFloat(data.full_mwt).toFixed(2)} g/mol`} />
        <Attribute name="Black box warning" value={data.bbw} />
        <Attribute name="Chirality" value={data.chirality} />
        <Attribute name="Max Phase" value={data.max_phase} tooltip="The maximum phase of development reached by this molecule." />
        <Attribute name="ALogP" value={data.alogp} />
        <Attribute name="Aromatic Rings" value={data.aromatic_rings} />
        <Attribute name="CX LogD pH7.4" value={data.cx_logd} />
        <Attribute name="CX LogP" value={data.cx_logp} />
        <Attribute name="CX Acidic pKa" value={data.cx_most_apka} />
        <Attribute name="Molecular Weight" value={data.full_mwt} />
        <Attribute name="Hydrogen Bond Acceptor" value={data.hba} />
        <Attribute name="Hydrogen Bond Donor" value={data.hbd} />
        <Attribute name="Heavy atoms" value={data.heavy_atoms} />
        <Attribute name="Polar Surface Area" value={data.psa} tooltip="The polar surface area (PSA) of a molecule is defined as the surface sum over all polar atoms or molecules, primarily oxygen and nitrogen, also including their attached hydrogen atoms." />
        <Attribute name="QED Weighted" value={data.qed_weighted} />
        <Attribute name="RO3 Pass" value={data.ro3_pass} />
        <Attribute name="#Rotatable Bonds" value={data.rtb} tooltip="The number of rotatable bonds (RBN) is the number of bonds which allow free rotation around themselves. These are defined as any single bond, not in a ring, bound to a nonterminal heavy atom." />
      </CardContent>
    </Card>
  );
};

export default ChEMBLCard;