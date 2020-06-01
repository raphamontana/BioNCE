import useSWR from 'swr';
import NGL from './NGLComponent';

const fetcher = url => fetch(url).then(r => r.json())

const MoleculeComponent = (props) => {
  const { data, error } = useSWR(`/api/smiles/${props.smiles}`, fetcher, { refreshInterval: 0 });
  if (error) return <div>failed to load</div>
  if (!data) return <div>loading...</div>

  let { pubchem } = data;
  console.log(data.bindings);

  return(
    <>
      <h2>Pubchem</h2>
      <p>CID: { pubchem.cid }</p>
      <p>SMILES: { data.smiles }</p>

      <h2>Binding data</h2>

      <hr/>


    </>
  );
};

export default MoleculeComponent;
//{data.bindings.map((structure, index) => (
//  <p key={structure}>{structure}</p>
//  ))}

//<span title="Download SDF" style="float: right;"><a href="https://files.rcsb.org/ligands/view/{{ data.pdb.id }}_model.sdf"><i className="ti-export"></i></a></span>