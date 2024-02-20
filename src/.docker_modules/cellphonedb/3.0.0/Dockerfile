FROM python:3.7-slim 

ENV RPY2_CFFI_MODE=ABI

RUN mkdir /app /docker
COPY requirements.txt /app
COPY requirements-web.txt /app
WORKDIR /app

RUN pip install -r requirements.txt
RUN pip install -r requirements-web.txt
RUN pip install cellphonedb==3.0.0


COPY . /app
COPY ./docker /docker
RUN mkdir -p /var/log/uwsgi
RUN touch /var/log/uwsgi/cellphonedb.log
RUN chmod +x /app

CMD ["/docker/run-system.sh"]
