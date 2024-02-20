FROM conda/miniconda2
MAINTAINER Laurent Modolo

ENV EMASE_VERSION=0.10.16
RUN conda init \
&& conda config --add channels r \
&& conda config --add channels bioconda \
&& conda create -y -n emase jupyter
SHELL ["conda", "run", "-n", "emase", "/bin/bash", "-c"]
RUN conda install -y -c kbchoi emase \
&& echo "conda activate emase" >> /root/.bashrc
RUN apt update && apt install -y procps

ENV PATH /usr/local/envs/emase/bin:/usr/local/condabin:$PATH
