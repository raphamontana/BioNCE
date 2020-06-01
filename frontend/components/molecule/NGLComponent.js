import { useEffect } from 'react';

function NGLComponent(props) {

  useEffect(() => {
    const loadNGL = async () => {
      const NGL = await import('ngl');
      const stage = new NGL.Stage("viewport");
      stage.loadFile("https://files.rcsb.org/download/4hhb.pdb", {defaultRepresentation: true});
    };
    loadNGL();
  },[]);

  return(
    <div id="viewport" style={{width: "400px", height:"300px"}} />
  );
}

export default NGLComponent;