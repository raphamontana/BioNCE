function NGLComponent({ viewportId, pdbFile }) {
  const [file, setFile] = React.useState(null);
  const [stage, setStage] = React.useState(null);

  if (pdbFile !== file) {
    setFile(pdbFile);
  }

  React.useEffect(() => {
    const loadNGL = async () => {
      const NGL = await import('ngl');
      const stage = new NGL.Stage(viewportId);
      stage.loadFile(file, {defaultRepresentation: true});
      setStage(stage);
    };
    loadNGL();
  },[]);

  React.useEffect(() => {
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