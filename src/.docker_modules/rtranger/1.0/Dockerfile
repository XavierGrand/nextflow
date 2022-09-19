FROM python:3.10.4-slim-buster

WORKDIR /app

COPY . . 

RUN apt update
RUN pip3 install --upgrade pip
RUN pip3 install -r requirements.txt

CMD ["python", "main.py", "--help"]