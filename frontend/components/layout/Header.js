import React from 'react';
import {
  AppBar, Avatar, Button, IconButton,
  InputBase, Link, Toolbar
} from '@material-ui/core';
import SearchIcon from '@material-ui/icons/Search';
import useStyles from './style';

const Header = () => {
  const classes = useStyles();
  const formRef = React.useRef();
  const [query, setQuery] = React.useState('');

  const submitOnEnter = (event) => {
    if (event.keyCode === 13) {
      formRef.current.submit();
    }
  };

  return (
    <AppBar position="absolute" className={classes.appBar}>
      <Toolbar paddingright="24">
        <Link href="/">
          <IconButton className={classes.menuButton} >
            <Avatar
              src={"/images/logomark.svg"}
              alt={"Bionce"}
            />
          </IconButton>
        </Link>
        <div className={classes.grow} />
        <div className={classes.search}>
          <form action="/molecules" method="POST" ref={formRef}>
            <input type="hidden" id="dataset" name="dataset" value={query} />
            <div className={classes.searchIcon}>
              <SearchIcon />
            </div>
            <InputBase
              value={query}
              onChange={e => setQuery(e.target.value)}
              placeholder="Search…"
              classes={{
                root: classes.inputRoot,
                input: classes.inputInput,
              }}
              inputProps={{ 'aria-label': 'search' }}
              onKeyDown={(e) => submitOnEnter(e)}
            />
          </form>
        </div>
        <div className={classes.grow} />
        <Button color="inherit" href="/help">
          Help
        </Button>
      </Toolbar>
    </AppBar>
  );
};

export default Header;