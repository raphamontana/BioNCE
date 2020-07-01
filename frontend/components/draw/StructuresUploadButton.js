import { DropzoneDialog } from 'material-ui-dropzone'
import Button from '@material-ui/core/Button';
import UploadIcon from '@material-ui/icons/CloudUpload';

const StructuresUploadButton = ({ addStructures }) => {
  const [open, setOpen] = React.useState(false);
  let fileReader;

  const handleFileRead = (e) => {
    const content = fileReader.result;
    addStructures(content);
  };

  const handleFileChosen = (files) => {
    setOpen(false);
    fileReader = new FileReader();
    fileReader.onloadend = handleFileRead;
    fileReader.readAsText(files[0]);
  };

  return (
    <>
      <Button
        variant="contained"
        color="primary"
        startIcon={<UploadIcon />}
        onClick={() => setOpen(true)}
      >
        Upload
      </Button>
      <DropzoneDialog
        open={open}
        onSave={handleFileChosen}
        acceptedFiles={['.csv', '.smi']}
        showPreviews={false}
        filesLimit={1}
        maxFileSize={1000000}
        onClose={() => setOpen(false)}
      />
    </>
  );
};

export default StructuresUploadButton;