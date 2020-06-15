-- Target is Cathepsin-L (CHEMBL3837) and Standard Type is IC50
SELECT m.chembl_id AS compound_chembl_id,
s.canonical_smiles,
r.compound_key,
COALESCE((d.pubmed_id::text), d.doi) AS pubmed_id_or_doi,
a.description                   AS assay_description,
act.standard_relation,
act.standard_value,
act.standard_units,
act.activity_comment,
t.chembl_id                    AS target_chembl_id,
t.pref_name                    AS target_name,
t.organism                     AS target_organism
FROM compound_structures s
  RIGHT JOIN molecule_dictionary m ON s.molregno = m.molregno
  JOIN compound_records r ON m.molregno = r.molregno
  JOIN docs d ON r.doc_id = d.doc_id
  JOIN activities act ON r.record_id = act.record_id
  JOIN assays a ON act.assay_id = a.assay_id
  JOIN target_dictionary t ON a.tid = t.tid
    AND t.chembl_id      = 'CHEMBL3837'
    AND act.standard_type = 'IC50'

    LIMIT 5;




-- Target is Cathepsin-L (CHEMBL3837) and standard_type is IC50
SELECT *
FROM compound_structures s
  RIGHT JOIN molecule_dictionary m ON s.molregno = m.molregno
  JOIN compound_records r ON m.molregno = r.molregno
  JOIN docs d ON r.doc_id = d.doc_id
  JOIN activities act ON r.record_id = act.record_id
  JOIN assays a ON act.assay_id = a.assay_id
  JOIN target_dictionary t ON a.tid = t.tid
    AND t.chembl_id      = 'CHEMBL3837'
    AND act.standard_type = 'IC50';

ALTER USER postgres WITH PASSWORD 'postgres';
