import { useState, useEffect } from 'react';

import parse from 'urlencoded-body-parser';
import Grid from '@material-ui/core/Grid';
import Layout from '../components/layout/Layout';
import MoleculeComponent from '../components/molecule/MoleculeComponent';

const Molecules = () => {

  const [molecules, setMolecules] = useState([]);

  useEffect(() => {
    let lstorage = window.localStorage.getItem('molecules');
    setMolecules(JSON.parse(lstorage));
    //TODO inserir um alerta de erro.
  },[]);

  return (
    <Layout>
      <Grid container spacing={3}>
        <Grid item xs>
          {molecules.map((option, index) => (
            <MoleculeComponent key={index} id={index} smiles={option} />
          ))}
        </Grid>
      </Grid>
    </Layout>
  );
};

export default Molecules;