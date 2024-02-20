FROM python:3.7
MAINTAINER Laurent Modolo

ENV KB_VERSION=0.24.4

RUN pip3 install kb_python==${KB_VERSION}
ENTRYPOINT ["/bin/sh"]
RUN ln -s /usr/local/lib/python3.7/site-packages/kb_python/bins/linux/bustools/bustools /usr/local/bin/
RUN ln -s /usr/local/lib/python3.7/site-packages/kb_python/bins/linux/kallisto/kallisto /usr/local/bin/
