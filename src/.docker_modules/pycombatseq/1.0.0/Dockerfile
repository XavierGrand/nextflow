FROM python:3.11-slim

# Installation des dépendances système nécessaires à la compilation de certains packages scientifiques
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Installation de inmoose et pandas
RUN pip install --no-cache-dir inmoose pandas

WORKDIR /data