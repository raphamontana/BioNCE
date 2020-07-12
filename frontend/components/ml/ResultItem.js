import useSWR from "swr";
import { LinearProgress, ListItem, ListItemIcon, ListItemText } from '@material-ui/core';
import { green, red } from '@material-ui/core/colors';
import CheckCircleIcon from '@material-ui/icons/CheckCircle';
import HighlightOffIcon from '@material-ui/icons/HighlightOff';

const fetcher = url => fetch(url).then(r => r.json());

const ResultItem = ({ id, model, smiles }) => {
  const { data, error } = useSWR(`/api/predict/${encodeURIComponent(model)}/${encodeURIComponent(smiles)}`, fetcher, { refreshInterval: 0 });
  if (error) return <div>Failed to predict: {smiles}</div>
  if (!data) return <LinearProgress />

  let icon;
  if (data.active === "Active") {
    icon = <CheckCircleIcon style={{ color: green[500] }} />;
  }
  else {
    icon = <HighlightOffIcon style={{ color: red.A400 }} />
  }
  return(
    <ListItem key={id}>
      <ListItemIcon>
        {icon}
      </ListItemIcon>
      <ListItemText
        primary={smiles}
        secondary={data.active}
      />
    </ListItem>
  );
}

export default ResultItem;