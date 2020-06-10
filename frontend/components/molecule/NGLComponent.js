import { useEffect, useState } from 'react';

function NGLComponent({ viewportId, pdbFile }) {
  const [file, setFile] = useState(null);
  const [stage, setStage] = useState(null);

  if (pdbFile !== file) {
    setFile(pdbFile);
  }

  useEffect(() => {
    const loadNGL = async () => {
      const NGL = await import('ngl');
      const stage = new NGL.Stage(viewportId);
      stage.loadFile(file, {defaultRepresentation: true});
      setStage(stage);
    };
    loadNGL();
  },[]);

  useEffect(() => {
    if (stage !== null) {
      stage.removeAllComponents();
      stage.loadFile(file, {defaultRepresentation: true});
    }
  },[file]);

  return(
    <div id={ viewportId } style={{ width: "100%", height:"250px" }} />
  );
}

export default NGLComponent;