FROM heroku/miniconda:3

# Grab requirements.txt.
ADD ./webapp/requirements.txt /tmp/requirements.txt

# Install dependencies
RUN pip install --upgrade pip
RUN pip install -qr /tmp/requirements.txt

# Add BioNCE code
ADD ./webapp /opt/webapp/
WORKDIR /opt/webapp
#ADD ./webapp/models /opt/webapp/models

RUN conda install -c conda-forge rdkit

CMD gunicorn --bind 0.0.0.0:$PORT wsgi