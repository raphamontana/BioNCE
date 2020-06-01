import Router from 'next/router';
import Button from '@material-ui/core/Button';
import Card from '@material-ui/core/Card';
import CardActions from '@material-ui/core/CardActions';
import CardContent from '@material-ui/core/CardContent';
import CardMedia from '@material-ui/core/CardMedia';
import Container from '@material-ui/core/Container';
import Typography from '@material-ui/core/Typography';
import Layout from '../components/layout/Layout';

const Custom404 = () => {
  return(
    <Layout>
      <Container maxWidth="sm">
        <Card>
          <CardMedia
            component="img"
            image='images/404.svg'
            height="100%"
            title="Something went wrong"
            alt="Something went wrong"
          />
          <CardContent>
            <Typography variant="body2" color="textSecondary" component="p">
              The page you are trying to reach does not exist or has been moved.
            </Typography>
          </CardContent>
          <CardActions>
            <Button size="small" color="primary" onClick={() => Router.back()}>
              Return
            </Button>
          </CardActions>
        </Card>
      </Container>
    </Layout>
  );
}

export default Custom404;