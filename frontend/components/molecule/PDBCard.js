import { useState } from 'react';
import { SimpleImg } from 'react-simple-img';
import { makeStyles } from '@material-ui/styles';
import { Avatar, Button, Card, CardActions, CardContent, CardHeader, FormControl,
         IconButton, InputLabel, MenuItem , Select, Tooltip } from '@material-ui/core';
import GetAppIcon from '@material-ui/icons/GetApp';
import LinkIcon from '@material-ui/icons/Link';
import Attribute from './Attribute';
import MolecularFormula from './MolecularFormula';

import NGL from './NGLComponent';

const useStyles = makeStyles((theme) => ({
  card: {
    minWidth: 290,
    maxWidth: 345,
    margin: theme.spacing(0.5),
  },
  formControl: {
    marginTop: theme.spacing(2),
    marginBottom: theme.spacing(2),
    width: "100%"
  },
}));

const PDBStructure = ({ structure }) => {
  const classes = useStyles();
  const [ligand, setLigand] = useState(structure.ligands[0]);

  return(
    <>
      <NGL
        pdbFile={ `rcsb://${ structure.id }` }
        viewportId={ "viewport" }
      />
      <p><b>Structure:</b> <a href={`https://www.rcsb.org/structure/${structure.id}`} target="_blank">{structure.id}</a></p>
      <p><b>Description:</b> {structure.description}</p>
      <p><b>Primary publication:</b>{' '}
        <a aria-label="Link to DOI" href={`http://doi.org/10.2210/pdb${structure.id}/pdb`} target="_blank">
          <LinkIcon />
        </a>
        <br />
        <b>Download:</b>{' '}
        <a aria-label="Link to download" href={`https://files.rcsb.org/download/${structure.id}.pdb`}>
          {structure.id}.pdb
        </a>
      </p>
      <FormControl variant="outlined" className={classes.formControl}>
        <InputLabel id="ligands-label">Ligands</InputLabel>
        <Select
          labelId="ligands-label"
          value={ligand}
          onChange={(e) => setLigand(e.target.value)}
          label="Ligands"
        >
          {structure.ligands.map((option, index) => (
            <MenuItem key={index} value={option} >{option}</MenuItem>
          ))}
        </Select>
      </FormControl>
      <Button
        size="small"
        color="primary"
        href={`https://www.rcsb.org/3d-view/${structure.id}?preset=ligandInteraction&sele=${ligand}`}
      >
        3D interaction
      </Button>
    </>
  );
}

const PDBCard = ({ data }) => {
  const classes = useStyles();
  const [structureIndex, setStructureIndex] = useState(0);

  return(
    <Card className={classes.card}>
      <CardHeader
        avatar={
          <a
            aria-label="PDB link"
            href={ `https://www.rcsb.org/ligand/${data.id}` }
            target="_blank"
          >
            <Avatar variant="rounded" alt="PDB icon" src="http://www.rcsb.org/favicon.ico" />
          </a>
        }
        action={
          <IconButton
            aria-label="Download SDF"
            href={`https://files.rcsb.org/ligands/view/${data.id}_model.sdf`}
          >
            <GetAppIcon />
          </IconButton>
        }
        title={ `Ligand ${data.id}` }
        subheader="Protein Data Bank"
      />
      <CardContent>
        <div align="center">
          <SimpleImg
            src={`https://cdn.rcsb.org/images/ccd/unlabeled/${data.id[0]}/${data.id}.svg`}
            alt="Compound structure"
            height={100}
          />
        </div>
        <Attribute name="Name" value={ data.name } />
        <MolecularFormula formula={ data.formula } />
        <Attribute name="Molecular Mass" value={ `${parseFloat(data.weight).toFixed(2)}` } />

        <FormControl variant="outlined" className={classes.formControl}>
          <InputLabel id="structures-label">Structures</InputLabel>
          <Select
            labelId="structures-label"
            value={structureIndex}
            onChange={(e) => setStructureIndex(e.target.value)}
            label="Structures"
          >
            {data.structures.map((option, index) => (
              <MenuItem key={index} value={index} >{option.id}</MenuItem>
            ))}
          </Select>
        </FormControl>
        <PDBStructure structure={data.structures[structureIndex]} />
      </CardContent>
    </Card>
  );
};

export default PDBCard;