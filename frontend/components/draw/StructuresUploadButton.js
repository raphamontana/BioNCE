import { useState, useEffect, useRef } from 'react'
import { DropzoneDialog } from 'material-ui-dropzone'
import Button from '@material-ui/core/Button';
import UploadIcon from '@material-ui/icons/CloudUpload';

const StructuresUploadButton = () => {
  const [open, setOpen] = useState(false);
  const [files, setFiles] = useState([]);
  const [smiles, setSmiles] = useState('');
  const formRef = useRef(null);

  const handleSave = async (files) => {
    setFiles(files);
    setOpen(false);
    const formData = new FormData();
    formData.append('smiFile', files[0]);
    try {
      let response = await fetch('/api/smiFileParser',
                                 { method: 'POST', body: formData });
      setSmiles(await response.json());
    }
    catch (err) {
      alert('Error:', err);
    };
  }

  useEffect(() => {
    if (smiles !== '') {
      formRef.current.submit();
    }
  }, [smiles]);

  return (
    <form action="/molecules" method="POST" ref={formRef}>
      <input type="hidden" id="structuresFile" name="structuresFile" value={smiles} />
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
        onSave={handleSave}
        acceptedFiles={['.csv', '.smi']}
        showPreviews={true}
        filesLimit={1}
        maxFileSize={500000}
        onClose={() => setOpen(false)}
      />
    </form>
  );
};

export default StructuresUploadButton;