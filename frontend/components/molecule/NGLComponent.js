import { useEffect } from 'react';

function NGLComponent({ viewportId, file }) {

  useEffect(() => {
    const loadNGL = async (viewportId) => {
      const NGL = await import('ngl');
      const stage = new NGL.Stage(viewportId);
      stage.loadFile(file, {defaultRepresentation: true});
    };
    loadNGL(viewportId);
  },[]);

  return(
    <div id={ viewportId } style={{ width: "100%", height:"250px" }} />
  );
}

export default NGLComponent;
//https://files.rcsb.org/download/4hhb.pdb