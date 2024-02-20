FROM ubuntu:20.04 

ENV DANPOS_VERSION="2f7f223"
ENV TZ=Europe/Paris
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone \
&& apt update -qq \
&& apt install --no-install-recommends software-properties-common -y dirmngr git python \
&& wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
&& add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN apt-get update \
&& apt-get install -y r-base=4.0.1 samtools \
&& git clone https://github.com/sklasfeld/DANPOS3.git \
&& cd DANPOS3 \
&& git checkout $DANPOS_VERSION \
&& pip install -r requirements.txt \
&& chmod +x DANPOS3/*.py

ENV PATH="/DANPOS3:${PATH}"
CMD [ "bash" ]
