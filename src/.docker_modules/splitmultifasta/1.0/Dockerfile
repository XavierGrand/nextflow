FROM python:3.10.4-buster

LABEL authors="Xavier Grand <xavier.grand@inserm.fr>"

RUN pip install --upgrade pip
RUN pip install biopython

COPY . /app

RUN useradd periodUser \
    && chown -R periodUser /app \
    && chmod +x /app/splitmultifasta.py
    
USER periodUser

CMD ["bash"]
