FROM alpine:3.12.0
MAINTAINER Lauret Modolo

ENV SABRE_VERSION=039a55e

ENV SABRE_DEPS g++ bash make zlib-dev git

RUN echo "http://ftp.acc.umu.se/mirror/alpinelinux.org/v3.11/main" > /etc/apk/repositories \
&& echo "http://ftp.acc.umu.se/mirror/alpinelinux.org/v3.11/community" >> /etc/apk/repositories \
&& sed -i -e 's/v[[:digit:]]\..*\//edge\//g' /etc/apk/repositories \
&& apk add --update --no-cache ${SABRE_DEPS} \
&& git clone https://github.com/najoshi/sabre.git \
&& cd sabre \
&& git checkout $SABRE_VERSION \
&& make \
&& mv sabre /usr/bin \
&& chmod +x /usr/bin/sabre


CMD ["bash"]
