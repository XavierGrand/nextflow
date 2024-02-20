FROM python:3.9-slim
MAINTAINER Lauret Modolo

ENV FLEXI_SPLITTER_VERSION=1.0.2
RUN apt-get update -qq \
    && apt-get install --no-install-recommends --yes \
        build-essential \
        procps
RUN pip3 install flexi-splitter==${FLEXI_SPLITTER_VERSION}
RUN apt-get remove --yes \
        build-essential

CMD [ "bash" ]