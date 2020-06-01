import AppBar from './AppBar';
import Sidebar from './Sidebar';
import useStyles from './style';
import Copyright from './Copyright';
import Container from '@material-ui/core/Container';

export default function Layout( { children } ) {
  const classes = useStyles();

  return (
    <div className={classes.root}>
      <AppBar />
      <Sidebar />
      <main
        flexgrow="1"
        height={"100vh"}
        overflow={"auto"}
        className={classes.content}>
        <div className={classes.appBarSpacer} />
        <Container maxWidth="lg" className={classes.container}>
          { children }
          <Copyright />
        </Container>
      </main>
    </div>
  );
}