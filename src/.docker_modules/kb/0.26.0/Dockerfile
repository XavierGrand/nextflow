FROM python:3.9-slim 

ENV KB_VERSION="0.26.0"

RUN apt update && apt install -y procps && pip3 install kb-python==${KB_VERSION}

COPY t2g.py /usr/bin/
COPY fix_t2g.py /usr/bin/

RUN chmod +x /usr/bin/t2g.py
RUN chmod +x /usr/bin/fix_t2g.py

CMD [ "bash" ]
