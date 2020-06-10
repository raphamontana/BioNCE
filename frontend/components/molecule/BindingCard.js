import ReactHtmlParser from 'react-html-parser';
import { makeStyles } from '@material-ui/styles';
import { Avatar, Card, CardContent, CardHeader, Grid, IconButton } from '@material-ui/core';
import GetAppIcon from '@material-ui/icons/GetApp';
import LinkIcon from '@material-ui/icons/Link';

const useStyles = makeStyles((theme) => ({
  card: {
    minWidth: 290,
    maxWidth: 345,
    margin: theme.spacing(0.5),
  }
}));

const InhibitionConstant = ({ measure, value, unit }) => {
  if (measure === "K<sub>a</sub>") {
    return null;
  }
  if (unit === "M^-1") {
    value = value / 10;
  }
  else if (unit === "mM") {
    value = value / 1000;
  }
  else if (unit === "uM") {
    value = value / 1000000;
  }
  else if (unit === "nM") {
    value = value / 1000000000;
  }
  let pValue = Math.log10(1 / value).toFixed(1);
  return(
    <>
      { ReactHtmlParser('p'+measure) } {" = " + pValue}
    </>
  );
}

const BindingInfo = ({ binding }) => {
  return(
    <Grid container justify="space-between">
      <Grid item>
        <b>{ binding.ligand }:</b> {"(Chain "} { binding.chain }{') '}
      </Grid>
      <Grid item>
        {ReactHtmlParser(binding.affinityMeasure)} { binding.relationType + ' ' + binding.affinityValue + ' ' + binding.affinityUnit }
        <br />
        <InhibitionConstant measure={binding.affinityMeasure} value={binding.affinityValue} unit={binding.affinityUnit} />
      </Grid>
    </Grid>
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
        title={ `Binding data for ${ligand}` }
        subheader="Binding MOAD"
      />
      <CardContent>
        {bindings.map((structure, index) => (
          <React.Fragment key={index}>
            <p>
              <b>Receptor:</b>{' '}
              <a href={`http://bindingmoad.org/pdbrecords/index/${ structure[0].proteinID }`} target="_blank">{ structure[0].proteinID }</a>
            </p>
            {structure.map((binding, index) => (
              <BindingInfo key={index} binding={binding} />
            ))}
            <hr />
          </React.Fragment>
        ))}
      </CardContent>
    </Card>
  );
}

export default BindingCard;