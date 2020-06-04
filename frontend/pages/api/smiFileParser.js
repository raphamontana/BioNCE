import { IncomingForm } from 'formidable';
import Papa from 'papaparse';
import fs from 'fs';

const SmiFileParser = async (req, res) => {
  const form = new IncomingForm();
  let parsed;
  try {
    const smiFile = await new Promise((resolve, reject) => {
      form.parse(req, (err, fields, files) => {
        if (err) return reject(err);
        if (typeof files.smiFile === 'undefined') {
          return reject("File not uploaded.");
        }
        resolve(files.smiFile.path);
      });
    });
    try {
      const smiFS = fs.readFileSync(smiFile, "utf8");
      parsed = Papa.parse(smiFS).data.map((value) => value[0]);
    }
    finally {
      fs.unlink(smiFile, (err) => {
        if (err) throw err;
      });
    }
  } catch (err) {
    return res.status(400).send(err.toString());
  }
  return res.status(200).json(parsed);
};

export const config = {
  api: {
    bodyParser: false,
  },
};

export default SmiFileParser;