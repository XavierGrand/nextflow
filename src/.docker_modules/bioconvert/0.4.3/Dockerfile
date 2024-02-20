FROM conda/miniconda3
MAINTAINER Laurent Modolo

ENV BIOCONVERT_VERSION="0.4.3"
RUN conda init \
&& conda config --add channels r \
&& conda config --add channels defaults \
&& conda config --add channels conda-forge \
&& conda config --add channels bioconda \
&& conda create -y -n bioconvert
SHELL ["conda", "run", "-n", "bioconvert", "/bin/bash", "-c"]
RUN conda install bioconvert \
&& echo "conda activate bioconvert" >> /root/.bashrc
RUN apt update && apt install -y procps

ENV PATH /usr/local/envs/bioconvert/bin:/usr/local/condabin:$PATH
