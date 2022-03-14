FROM heroku/miniconda:3

# Grab requirements.txt.
ADD ./webapp/requirements.txt /tmp/requirements.txt

# Install dependencies
RUN pip install --upgrade pip
RUN pip install -qr /tmp/requirements.txt

# Add BioNCE code
ADD ./webapp /opt/webapp/
WORKDIR /opt/webapp
#ADD ./webapp/endpoints /opt/webapp/endpoints
#ADD ./webapp/models /opt/webapp/models
#ADD ./webapp/ad /opt/webapp/ad
#ADD ./webapp/pains /opt/webapp/pains


RUN conda install -c rdkit rdkit
RUN conda install -c chembl chembl_structure_pipeline

CMD gunicorn --bind 0.0.0.0:$PORT wsgi
