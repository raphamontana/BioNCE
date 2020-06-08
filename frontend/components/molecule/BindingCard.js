import ReactHtmlParser from 'react-html-parser';
import { makeStyles } from '@material-ui/styles';
import { Avatar, Card, CardContent, CardHeader, FormControl,
         IconButton, InputLabel, MenuItem , Select, Tooltip } from '@material-ui/core';
import GetAppIcon from '@material-ui/icons/GetApp';
import LinkIcon from '@material-ui/icons/Link';

import NGL from './NGLComponent';

const useStyles = makeStyles((theme) => ({
  card: {
    minWidth: 290,
    maxWidth: 345,
    margin: theme.spacing(0.5),
  }
}));

const BindingInfo = ({ binding }) => {
  return(
    <p>
      <b>{ binding.ligand }</b>
      {"(Chain "} { binding.chain }{') '}
      <span> {ReactHtmlParser(binding.affinityMeasure)} { binding.relationType + ' ' + binding.affinityValue + ' ' + binding.affinityUnit }</span>
    </p>
  );
}

const BindingCard = ({ ligand, bindings }) => {
  const classes = useStyles();

  return(
    <Card className={classes.card}>
      <CardHeader
        avatar={
          <a
            aria-label="Binding link"
            href={ `http://bindingmoad.org/Search/showsearch/%2A/%2A/%2A/%2A/${ligand}/%2A` }
            target="_blank"
          >
            <Avatar variant="rounded" alt="Binding icon" src="http://bindingmoad.org/img/blockm.png" />
          </a>
        }
        action={
          <IconButton
            aria-label="Download SDF"
          >
            <LinkIcon />
          </IconButton>
        }
        title={ `Binding data ${ligand}` }
        subheader="Binding MOAD"
      />
      <CardContent>
        <NGL
          file={ `https://files.rcsb.org/download/${ bindings[0][0].proteinID }.pdb` }
          viewportId={ "viewport-" + ligand }
        />
        {bindings.map((structure, index) => (
          <React.Fragment key={index}>
            <p>
              <b>Receptor:</b>{' '}
              <a href={`http://bindingmoad.org/pdbrecords/index/${ structure[0].proteinID }`} target="_blank">{ structure[0].proteinID }</a>
            </p>
            {structure.map((binding, index) => (
              <BindingInfo key={index} binding={binding} />
            ))}
          </React.Fragment>
        ))}
      </CardContent>
    </Card>
  );
}

export default BindingCard;