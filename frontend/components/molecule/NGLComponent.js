import { useEffect } from 'react';

function NGLComponent(props) {

  useEffect(() => {
    const loadNGL = async (viewportId) => {
      const NGL = await import('ngl');
      const stage = new NGL.Stage(viewportId);
      stage.loadFile("https://files.rcsb.org/download/4hhb.pdb", {defaultRepresentation: true});
    };
    loadNGL(props.viewportId);
  },[]);

  return(
    <div id={ props.viewportId } style={{ width: "400px", height:"300px" }} />
  );
}

export default NGLComponent;