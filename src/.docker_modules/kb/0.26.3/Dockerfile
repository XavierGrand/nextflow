FROM python:3.9-slim 

ENV KB_VERSION="0.26.3"

RUN apt update && apt install -y procps make gcc zlib1g-dev libbz2-dev  	libcurl4 liblzma-dev \
        && pip3 install pysam anndata h5py Jinja2 loompy nbconvert nbformat ngs-tools numpy pandas plotly scanpy scikit-learn tqdm \
        && pip3 install kb-python==${KB_VERSION} gffutils

COPY t2g.py /usr/bin/

RUN chmod +x /usr/bin/t2g.py

CMD [ "bash" ]
