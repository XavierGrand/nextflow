FROM jcrist/alpine-conda:4.6.8

RUN /opt/conda/bin/conda config --add channels r \
    && /opt/conda/bin/conda config --add channels bioconda \
    && /opt/conda/bin/conda install --yes \
        -c kbchoi g2gtools pycparser setuptools\
        nomkl \
    && /opt/conda/bin/conda clean -afy \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete

